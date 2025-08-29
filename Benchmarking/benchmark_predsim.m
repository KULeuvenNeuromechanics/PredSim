function [] = benchmark_predsim(S,osim_path,S_benchmark)
%benchmark_predsim function to benchmark your predsim simulation workflow
% as in afschrift 2025 (benchmarking the predictive capability of human gait
% simulations)
%   Input arguments:
%       (1) S: this is the matlab structure with settings that you use in
%       your predsim simulation workfloz
%       (2) osim_path: path to the opensim model (same in in predsim
%       workflow)
%       (3) S_benchmark: specific settings for benchmarking (see readme)


% also keep diary of this script ? does it work to have two diary function
% inside each other ?
log_name = fullfile(S_benchmark.out_folder,'log_benchmark.txt');
diary(log_name);

% get casadi path
if ~isfield(S.solver,'CasADi_path')
    try
        S.solver.CasADi_path = casadi.GlobalOptions.getCasadiPath();
    catch
        error("Please add CasADi to the matlab search path, or pass the path " + ...
            "to your CasADi installation (top folder) to S.solver.CasADi_path.")
    end
elseif ~isempty(S.solver.CasADi_path) && ~isfolder(S.solver.CasADi_path)
    error("Unable to find the path assigned to S.solver.CasADi_path:" + ...
        " \n\t%s",S.solver.CasADi_path)
end


% The only input is S and osim_path, maybe additional input S_benchmark
S_input = S; % make copy of S
osim_path_default = osim_path;
if ~isfolder(S_benchmark.out_folder)
    mkdir(S_benchmark.out_folder)
end
S.misc.save_folder = S_benchmark.out_folder;


% get path information

% flow of the function
%-----------------------
% 0. pre-processing default model
%       I think the most clean way is to add additional settings to control
%       the flow in the run_pred_sim script
S.flow_control.pre_processing_only = true;
S.solver.run_as_batch_job = false;
run_pred_sim(S,osim_path_default);
diary(log_name);

%
% 1. check if we have to make multiple .dll files
%    (this is needed for walking on a slope and with added mass)
%    if this is the case:
%       - make new .osim files
%       - create new .dll files
%       - copy preprocessing information from default model to new models
%           - F_Falisse_et_al_2022_IO.mat
%           - f_lMT_vMT_dM_poly_3_9
%
studies_convertmodel_needed = {'koelewijn2019','browning2008',...
    'schertzer2014','abe2015','gomenuka2013','mcdonald2022','huang2014'};
% check if we have to convert a model
bool_convertmodels = false;
for istudy = 1:length(S_benchmark.studies)
    if any(strcmp(S_benchmark.studies{istudy},studies_convertmodel_needed))
        bool_convertmodels = true;
    end
end
% convert models if needed
if bool_convertmodels
    % message
    disp(' ')
    disp('started converting models for benchmarking studies');
    disp(' ')

    % store all subjects in benchmarking folder
    S_benchmark.converted_models = [];

    %       koelewijn2019
    if any(strcmp(S_benchmark.studies,'koelewijn2019'))
        % desired models
        ct = 1; % counter models
        slopes = [-8 8];
        for islope = 1:length(slopes)
            % create osim model with adapted gravity
            if slopes(islope)<0
                model_name = [S.subject.name '_down_' num2str(abs(slopes(islope)))];
            elseif slopes(islope)>0
                model_name = [S.subject.name '_up_' num2str(abs(slopes(islope)))];
            end
            out_folder =  fullfile(S.misc.main_path,'Subjects',model_name);
            out_modelname = fullfile(out_folder,[model_name '.osim']);
            adapt_gravity_model(osim_path_default,slopes(islope),out_modelname);
            S_benchmark.converted_models.koelewijn2019.modelnames{ct} = model_name;
            S_benchmark.converted_models.koelewijn2019.osim_path{ct} = out_modelname;

            % copy model geometry file and settingsfile
            copy_musclegeom_information(osim_path_default,out_modelname);
            copy_modelsettingsfile(osim_path_default,out_modelname)

            % convert model
            S_temp = S_input;
            S_temp.flow_control.pre_processing_only = true;
            S_temp.solver.run_as_batch_job = false;
            S_temp.subject.name = model_name;
            S_temp.misc.save_folder = S_benchmark.out_folder;
            run_pred_sim(S_temp,out_modelname);
            diary(log_name);
            ct = ct + 1;
        end
        % add default model (level walking)
        S_benchmark.converted_models.koelewijn2019.modelnames{ct}= S_input.subject.name;
        S_benchmark.converted_models.koelewijn2019.osim_path{ct}= osim_path_default;
    end

    % Browning
    if any(strcmp(S_benchmark.studies,'browning2008'))
        % create osim model for browning study
        [browning2008] = adapt_model_Browning(S,osim_path_default);
        S_benchmark.converted_models.browning2008 = browning2008;
        % add default model (walking without added mass)
        ct = length(S_benchmark.converted_models.browning2008.modelnames) +1;
        S_benchmark.converted_models.browning2008.modelnames{ct}= S_input.subject.name;
        S_benchmark.converted_models.browning2008.osim_path{ct}= osim_path_default;
        % create dlls
        for imodel  = 1:length(browning2008.osim_path)
            % copy the muscle-tendon information
            copy_musclegeom_information(osim_path_default,browning2008.osim_path{imodel});
            copy_modelsettingsfile(osim_path_default,browning2008.osim_path{imodel})

            % convert model
            S_temp = S_input;
            S_temp.flow_control.pre_processing_only = true;
            S_temp.solver.run_as_batch_job = false;
            S_temp.subject.name = browning2008.modelnames{imodel};
            S_temp.misc.save_folder = S_benchmark.out_folder;
            run_pred_sim(S_temp,browning2008.osim_path{imodel});
            diary(log_name);
        end
        
    end
end




% 2. Run batch simulations

% add check if simulation result already exists, if this is the case do not
% run the simulation
disp('')
disp('preprocessing of all models finished')
disp(' ')
disp('start simulations');
disp(' ')


% batch replication of specific studies
if isfield(S_benchmark,'studies') && ~isempty(S_benchmark.studies)

    if any(strcmp(S_benchmark.studies,'vanderzee2022'))
        % ----- Van Der Zee 2022 -------
        % DOI: 10.1038/s41597-022-01817-1
        % we currently only benchmark the gait speeds without imposed step
        % frequency or imposed step width
        % these are the trials with speeds and id:
        gait_speeds = [0.7 0.9 1.1 1.6 1.8 2 1.4];
        id_trials = [2 5 8 11 14 16 32]; % see publication
        S_benchmark.vanderzee.gait_speeds = gait_speeds;
        S_benchmark.vanderzee.id_trials = id_trials;
        S_benchmark.vanderzee.names = cell(length(id_trials),1);
        S_benchmark.vanderzee.slopes = 0;
        S_benchmark.vanderzee.addedmass = 0;
        for i_speed = 1:length(gait_speeds)
            % start from default input settings
            S = S_input;
            % set forward velocity
            S.misc.forward_velocity = gait_speeds(i_speed);
            % adapt folder to save results
            S.misc.save_folder  = fullfile(S_benchmark.out_folder,'vanderzee2022',...
                ['speed_' num2str(round(gait_speeds(i_speed)*100))...
                '_id_' num2str(id_trials(i_speed))]);
            % check if save_folder already exists and contains a matfile,
            % if this is the case do not run the simulation
            mat_files = dir(fullfile(S.misc.save_folder,'*.mat'));
            if isempty(mat_files)
                % run predsim
                runPredSim(S, osim_path_default);
            end
            % append name
            S_benchmark.vanderzee.names{i_speed} = ['speed_' ...
                num2str(round(gait_speeds(i_speed)*100)) '_id_' ...
                num2str(id_trials(i_speed))];
        end
    end

    % run simulations as in experiment koelewijn2019
    if any(strcmp(S_benchmark.studies,'koelewijn2019'))
        % specify all conditions
        S_benchmark.koelewijn.gait_speeds = [0.8 1.3];
        S_benchmark.koelewijn.names = cell(length(id_trials),1);
        S_benchmark.koelewijn.slopes = [-8 0 8];
        S_benchmark.koelewijn.addedmass = 0;
        % the idea is always to run for all models associated with the
        % study the range of gait speeds
        ct_sim = 1;
        for imodel = length(S_benchmark.converted_models.koelewijn2019.modelnames)
            % name of the model (model on slope ?
            model_name = S_benchmark.converted_models.koelewijn2019.modelnames{imodel};
            osim_path_sel = S_benchmark.converted_models.koelewijn2019.osim_path{imodel};
            for i_speed = 1:length(S_benchmark.koelewijn.gait_speeds)
                % copy input settings
                S = S_input;
                % set forward velocity
                S.misc.forward_velocity = S_benchmark.koelewijn.gait_speeds(i_speed);
                % adapt folder to save results
                save_name  =[model_name, '_speed_' num2str(round(gait_speeds(i_speed)*100))];
                S.misc.save_folder  = fullfile(S_benchmark.out_folder,'koelewijn2019',...
                    save_name);
                % adapt subject
                S.subject.name = model_name;
                % check if save_folder already exists and contains a matfile,
                % if this is the case do not run the simulation
                mat_files = dir(fullfile(S.misc.save_folder,'*.mat'));
                if isempty(mat_files)                
                    % run predsim
                    runPredSim(S, osim_path_sel);
                end
                % append name
                S_benchmark.koelewijn.names{ct_sim} = save_name;
                ct_sim = ct_sim+1;
            end
        end
    end

    % run simulations as in browning 2008
    if any(strcmp(S_benchmark.studies,'browning2008'))        
        S_benchmark.browning.gait_speeds = [0.8 1.3];
        S_benchmark.browning.names = cell(length(id_trials),1);
        S_benchmark.browning.slopes = 0;
        S_benchmark.browning.addedmass = 0;
        ct_sim = 1;
        for imodel = length(S_benchmark.converted_models.browning2008.modelnames)
            % name of the model (model on slope ?
            model_name = S_benchmark.converted_models.browning2008.modelnames{imodel};
            osim_path_sel = S_benchmark.converted_models.browning2008.osim_path{imodel};
            % copy input settings
            S = S_input;
            % set forward velocity
            S.misc.forward_velocity = 1.25;
            % adapt folder to save results
            save_name  = model_name;
            S.misc.save_folder  = fullfile(S_benchmark.out_folder,'browning2008',...
                save_name);
            % adapt subject
            S.subject.name = model_name;
            % check if save_folder already exists and contains a matfile,
            % if this is the case do not run the simulation
            mat_files = dir(fullfile(S.misc.save_folder,'*.mat'));
            if isempty(mat_files)
                % run predsim
                runPredSim(S, osim_path_sel);
            end
            % append name
            S_benchmark.browning2008.names{ct_sim} = save_name;
            ct_sim = ct_sim+1;

        end


    end
end

% batch simulation gait speeds
if isfield(S_benchmark,'gait_speeds') && S_benchmark.gait_speeds

    % get set of gait speeds to test
    if isfield(S_benchmark,'gait_speeds_selection') && ...
            ~isempty(S_benchmark.gait_speeds_selection)
        gait_speeds = S_benchmark.gait_speeds_selection;
    else
        if isfield(S_benchmark,'gait_speed_range') && ...
                ~isempty(S_benchmark.gait_speed_range)
            speed_min = min(S_benchmark.gait_speed_range);
            speed_max = max(S_benchmark.gait_speed_range);
            gait_speeds = speed_min:0.1:speed_max;
        end
    end
    if isempty(gait_speeds)
        gait_speeds = 0.6:0.1:2;
    end

    % run predictive simulations
    S_benchmark.gaitspeed.gait_speeds = gait_speeds;
    S_benchmark.gaitspeed.names = cell(length(gait_speeds),1);
    for i_speed = 1:length(gait_speeds)
        % start from default input settings
        S = S_input;
        % set forward velocity
        S.misc.forward_velocity = gait_speeds(i_speed);
        % adapt folder to save results
        S.misc.save_folder  = fullfile(S_benchmark.out_folder,'gait_speeds',...
            ['speed_' num2str(round(gait_speeds(i_speed)*100))]);
        S_benchmark.gaitspeed.names{i_speed} = ...
            ['speed_' num2str(round(gait_speeds(i_speed)*100))];
        % check if save_folder already exists and contains a matfile,
        % if this is the case do not run the simulation
        mat_files = dir(fullfile(S.misc.save_folder,'*.mat'));
        if isempty(mat_files)
            % run predsim
            runPredSim(S, osim_path_default);
        end

    end
end

% save all settings for the benchmark procedure in the output folder. We
% will use this settings once all simulations are finished to benchmark the
% simulation results / compare it to experimental data
S = S_input;
save(fullfile(S_benchmark.out_folder,'benchmark_settings.mat'),'S',...
    'S_benchmark','osim_path');


end