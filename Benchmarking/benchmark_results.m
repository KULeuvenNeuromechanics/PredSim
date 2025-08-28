function [] = benchmark_results(benchmarking_folder,varargin)
%benchmark_results Compares simulations to experiments
%   input arguments:
%       (1) benchmarking_folder = folder with benchmarking results. this is
%       the S_benchmark.out_folder folder when running benchmark_predsim
%       (2) optional input arguments:
%           - 'BoolPlot', true, makes some default plots with benchmarking results


% Note:
% this is work in progress. the plan is to work out a minimal example with
% benchmarking the study of Van Der Zee with variations in gait speed. I
% will discuss this example at a PredSim meeting


%% input parser
p = inputParser;
addParameter(p, 'BoolPlot', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'dofs', {}, @(x) iscellstr(x) || isstring(x));
parse(p, varargin{:});
BoolPlot = logical(p.Results.BoolPlot);
dofs_plot = cellstr(p.Results.dofs);

% defaults settings
if isempty(dofs_plot)
    dofs_plot = {'ankle_angle_r','knee_angle_r','hip_flexion_r'};
end
nsubplot_dofs = length(dofs_plot);


%% load benchmarking settings

load(fullfile(benchmarking_folder,'benchmark_settings.mat'),'S','osim_path','S_benchmark');


%% load experimental data

% DOI 10.5281/zenodo.14524619
data_exp_path = fullfile(benchmarking_folder,'exp_data');
if ~isfolder(data_exp_path)
    mkdir(data_exp_path);
end

% general data
% check if we have to download the data
if ~isfolder(fullfile(benchmarking_folder,'exp_data','general_exp_data'))
    % download data as in Afschrift et al. 2025
    LinkZipFile = ['https://zenodo.org/records/14524620/files/',...
        'PredSim_gait_conditions-master_19_12_2024.zip?download=1'];
    websave(fullfile(data_exp_path,'githubrepo_data.zip'),LinkZipFile);
    unzip(fullfile(data_exp_path,'githubrepo_data.zip'),...
        fullfile(data_exp_path,'github_repo_temp'));
    delete(fullfile(data_exp_path,'githubrepo_data.zip'));
    %       change data format a bit here, we don't need everything
    copyfile(fullfile(data_exp_path,'github_repo_temp',...
        'PredSim_gait_conditions-master','ExperimentalData'),...
        fullfile(data_exp_path,'general_exp_data'))
    rmdir(fullfile(data_exp_path,'github_repo_temp'),'s');
end

% Van Der Zee kinematic data
if isfield(S_benchmark,'studies') && ~isempty(S_benchmark.studies)
    if any(strcmp(S_benchmark.studies,'vanderzee2022'))
        if ~isfolder(fullfile(benchmarking_folder,'exp_data','vanderzee2022'))
            % url to datafile
            LinkZipFile = ['https://zenodo.org/records/14524620/files/',...
                'VanDerZee_Averages_Addbiomech.zip?download=1'];
            % download data
            websave(fullfile(data_exp_path,'vanderzee2022.zip'),LinkZipFile);
            % unzip and delete zip file
            unzip(fullfile(data_exp_path,'vanderzee2022.zip'),...
                fullfile(data_exp_path,'vanderzee2022'));
            delete(fullfile(data_exp_path,'vanderzee2022.zip'));
        end
    end
end

%% Attach experimental data to simulation results

if isfield(S_benchmark,'studies') && ~isempty(S_benchmark.studies)

    % subject properties
    Studyname = {'Gomenuka','Browning','Schertzer','Huang','Koelewijn',...
        'vanderzee2022','McDonald','Abe2015','Strutzenberger'};
    prop_leg_length = 0.55;
    Mass = [71.6, 74.16, 74.88, 71.1, 70, 73.5, 69.6, 58.9, 73.1];
    Height = [1.78, 1.82, 1.78, 0.99/prop_leg_length, 1.73, 1.76, 1.70, 1.70, 1.77];
    LegLength = Height.*prop_leg_length;
    Exp_SubjProp = table(Studyname',Mass',Height',LegLength','VariableNames',{'Study','mass','height','LegLength'});
    msim = 62;
    Height_Sim = 1.70;
    Lsim = Height_Sim.*prop_leg_length;
    g = 9.81;



    if any(strcmp(S_benchmark.studies,'vanderzee2022'))
        % load the sf datafile, this contains unfortunatly all conditions
        % todo: see in script of benchmarking paper to select correct ones
        sf_datafile = fullfile(data_exp_path,'vanderzee2022','stridefreq.csv');
        % loop over simulations
        n_sim = length(S_benchmark.vanderzee.names);
        for isim = 1:n_sim
            % path to experimental data
            speed = S_benchmark.vanderzee.gait_speeds(isim);
            grf_datafile = fullfile(data_exp_path,'vanderzee2022',...
                ['mean_' num2str(speed*10) '_GRF.csv']);
            grf_datafile_std = fullfile(data_exp_path,'vanderzee2022',...
                ['mean_' num2str(speed*10) '_GRF_std.csv']);
            id_datafile = fullfile(data_exp_path,'vanderzee2022',...
                ['mean_' num2str(speed*10) '_ID.csv']);
            id_datafile_std = fullfile(data_exp_path,'vanderzee2022',...
                ['mean_' num2str(speed*10) '_ID_std.csv']);
            ik_datafile = fullfile(data_exp_path,'vanderzee2022',...
                ['mean_' num2str(speed*10) '_IK.csv']);
            ik_datafile_std = fullfile(data_exp_path,'vanderzee2022',...
                ['mean_' num2str(speed*10) '_IK_std.csv']);


            % attach to simulation results
            % store everything in a matlab datastruture with the fields
            %   - grf           OR: GRF_r and GRF_l
            %   - grf_std       OR: GRF_r_std and GRF_l_std
            %   -
            %
            %   - ik
            %   - ik_std
            %   - id
            %   - id_std
            %   - stride_frequency
            %   - Pmetab_mean
            %   - subject_height: (average) height of participant
            %   - subject_mass: (average) mass of participant
            % note that for ik and id it is important that the same model is used
            % for experiment and predictive simulation. The generalized coordinates
            % should at least have the same name between both models if you want that
            % the code runs.
            % I decided to have everything here with units (so not non-dim or normalised)
            % lets hope I don't regret this decision.

            sim_res_folder = fullfile(benchmarking_folder,'vanderzee2022',...
                S_benchmark.vanderzee.names{isim});
            mat_files = dir(fullfile(sim_res_folder,'*.mat'));
            if length(mat_files) ~= 1
                disp('warning mutiple mat files in folder')
                disp(sim_res_folder);
                disp(['assumes that file ' mat_files(1).name, ...
                    'contains the simulation results'])
            end
            sim_res_file = fullfile(mat_files(1).folder, mat_files(1).name);

            benchmark.grf = readtable(grf_datafile); % check if already norm to mass
            benchmark.grf_std = readtable(grf_datafile_std); % check if already norm to mass
            benchmark.id = readtable(id_datafile); % check if already norm to mass
            benchmark.id_std = readtable(id_datafile_std);
            benchmark.ik= readtable(ik_datafile);
            benchmark.ik_std = readtable(ik_datafile_std);
            benchmark.Pmetab_mean = NaN; % no data
            Lexp = Exp_SubjProp.height(strcmp(Exp_SubjProp.Study,'vanderzee2022'));
            mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'vanderzee2022'));
            benchmark.subject_height = Lexp;
            benchmark.subject_mass = mexp; % to do check in publication

            % get average stride frequency across subjects (not yet done in
            % datafile)

            sf_all = readtable(sf_datafile);
            sf_mean = nanmean(sf_all.stridefreq(sf_all.walkspeed == speed));
            sf_mean = sf_mean*(Lexp*0.5); % neede because we exported norm stride frequency with a weird norm for some reason
            benchmark.stride_frequency = sf_mean;

            % adapt some data as grf and moments are in N/kg
            % adapt all cols expect the time col
            id_names = benchmark.id.Properties.VariableNames;
            for idof = 1:length(id_names)
                if ~any(strcmp(id_names{idof},{'time','Time'}))
                    benchmark.id.(id_names{idof}) = ...
                        benchmark.id.(id_names{idof})*benchmark.subject_mass;
                    benchmark.id_std.(id_names{idof}) = ...
                        benchmark.id_std.(id_names{idof})*benchmark.subject_mass;
                end
            end

            grf_names = benchmark.grf.Properties.VariableNames;
            for igrf = 1:length(grf_names)
                if ~any(strcmp(grf_names{igrf},{'time','Time'}))
                    benchmark.grf.(grf_names{igrf}) = ...
                        benchmark.grf.(grf_names{igrf}) * benchmark.subject_mass;
                end
            end

            % note I decide here on grf headers like this:
            %  Flx         Fly        Flz         Frx       Fry        Frz
            % I might want to do this with the same headers as in predsim
            % grf_headers = GRF_r: [100×3 double] and GRF_l: [100×3 double]
            % ToDo: discuss this in predsim meeting. I think I'm in favour of adapting
            % experimental data such that it has the same structure as the predsim results
            benchmark.GRF_r = [-benchmark.grf.Flx, benchmark.grf.Fly, -benchmark.grf.Flz]; % have to switch l and r and rotate here for some reason
            benchmark.GRF_l = [-benchmark.grf.Frx, benchmark.grf.Fry, -benchmark.grf.Frz]; % have to switch l and r and rotate here for some reason
            benchmark.GRF_r_std = [benchmark.grf_std.Frx, benchmark.grf_std.Fry, benchmark.grf_std.Frz];
            benchmark.GRF_l_std = [benchmark.grf_std.Flx, benchmark.grf_std.Fly, benchmark.grf_std.Flz];

            % add benchmark to .mat file
            save(sim_res_file, 'benchmark', "-append")

        end
    end
end


%% Default plots
if BoolPlot
    % allowed to adapt figure defaults ?
    set(0,'defaultLineLineWidth',1.6);
    set(0,'defaultAxesLineWidth',1.2);
    set(0,'defaultAxesFontSize',12);
    set(0,'defaultFigureColor',[1 1 1]);
    set(0,'defaultLineMarkerSize',4);
    set(0,'defaultAxesBox','off');

    if any(strcmp(S_benchmark.studies,'vanderzee2022'))
        % default function to plot van der zee results
        % to do: add this to a function
        %
        speeds = S_benchmark.vanderzee.gait_speeds;        
        [speeds_sort,isort] = sort(speeds);
        for i= 1:length(speeds_sort)
            WalkSpeed_legend{i} = [num2str(speeds_sort(i)) 'ms^{-1}'];
        end
        % Colors_Speeds = [146,133,101;... % 0.7
        %     126,108,62;... % 0.9
        %     238,202,102;...% 1.1
        %     125,162,197;...% 1.6
        %     238,202,102;...% 1.8
        %     125,162,197;...% 2.0            
        %     238,202,102]./255 % 1.4; % these colors are pretty bad, also in paper ?
        % Colors_Speeds = flipud(Colors_Speeds);
        Colors_Speeds = copper(length(speeds)+2);
        % load data of all simulation conditions first
        for ispeed_sort = 1:length(speeds_sort)
            ispeed = isort(ispeed_sort);
            sim_res_folder = fullfile(benchmarking_folder,'vanderzee2022',...
                S_benchmark.vanderzee.names{ispeed});
            mat_files = dir(fullfile(sim_res_folder,'*.mat'));
            if length(mat_files) ~= 1
                disp('warning mutiple mat files in folder')
                disp(sim_res_folder);
                disp(['assumes that file ' mat_files(1).name, ...
                    'contains the simulation results'])
            end
            sim_res_file = fullfile(mat_files(1).folder, mat_files(1).name);
            load(sim_res_file,'R','benchmark');
            Dat(ispeed_sort).benchmark = benchmark;
            Dat(ispeed_sort).R = R;
        end

        % scale factors to nondim everything
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'vanderzee2022'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'vanderzee2022'));
        % Step frequency was normalized by \sqrt{g/l}, 
        % joint moments by mgl, ground reaction forces by mg, 
        % and (metabolic) power by mg^{1.5}l^{0.5}.

        % plot kinematics
        figure('Name','vanderzee: kinematics','Color',[1 1 1]);
        
        for ispeed = 1: length(speeds_sort)
            for idof = 1:length(dofs_plot)
                % plot experimental data
                subplot(2, nsubplot_dofs, idof)
                plot(Dat(ispeed).benchmark.ik.(dofs_plot{idof})*180/pi,...
                    'Color',Colors_Speeds(ispeed,:)); hold on;


                % plot simulation data
                subplot(2, nsubplot_dofs, idof+nsubplot_dofs)
                icol = strcmp(dofs_plot{idof},Dat(ispeed).R.colheaders.coordinates);
                dsel = Dat(ispeed).R.kinematics.Qs(:,icol);
                dsel_int = interp1(1:length(dsel),dsel,linspace(1,length(dsel),100));
                l(ispeed) = plot(dsel_int,'Color',Colors_Speeds(ispeed,:)); hold on;
            end
        end
        legend(l,WalkSpeed_legend,'Interpreter','tex');
        legend boxoff;
        for isubpl =1:length(dofs_plot)*2
            subplot(2, nsubplot_dofs,isubpl)
            set(gca,'box','off');
            set(gca,'FontSize',10);
            if isubpl<=length(dofs_plot)
                title(dofs_plot{isubpl},'interpreter','none');
            else
                xlabel('% gait cycle');
            end
            if isubpl == 1
                ylabel({'experiment','joint angle [deg]'});
            elseif isubpl == (length(dofs_plot)+1)
                ylabel({'simulation','joint angle [deg]'});
            end
        end

        % plot joint moments
        figure('Name','vanderzee: kinetics','Color',[1 1 1]);
        
        for ispeed = 1: length(speeds_sort)
            for idof = 1:length(dofs_plot)
                % plot experimental data
                subplot(2, nsubplot_dofs, idof)
                id_exp = Dat(ispeed).benchmark.id.(dofs_plot{idof});
                id_exp = id_exp./(mexp*g*Lexp).*(msim*g*Lsim); % scale to subject
                plot(id_exp,'Color',Colors_Speeds(ispeed,:)); hold on;


                % plot simulation data
                subplot(2, nsubplot_dofs, idof+nsubplot_dofs)
                icol = strcmp(dofs_plot{idof},Dat(ispeed).R.colheaders.coordinates);
                dsel = Dat(ispeed).R.kinetics.T_ID(:,icol);
                dsel_int = interp1(1:length(dsel),dsel,linspace(1,length(dsel),100));
                l(ispeed) = plot(dsel_int,'Color',Colors_Speeds(ispeed,:)); hold on;
            end
        end
        legend(l,WalkSpeed_legend,'Interpreter','tex');
        legend boxoff;
        for isubpl =1:length(dofs_plot)*2
            subplot(2, nsubplot_dofs,isubpl)
            set(gca,'box','off');
            set(gca,'FontSize',10);
            if isubpl<=length(dofs_plot)
                title(dofs_plot{isubpl},'interpreter','none');
            else
                xlabel('% gait cycle');
            end
            if isubpl == 1
                ylabel({'experiment','joint moment [Nm]'});
            elseif isubpl == (length(dofs_plot)+1)
                ylabel({'simulation','joint moment [Nm]'});
            end
        end


        % plot ground reaction forces
        figure('Name','vanderzee: grf','Color',[1 1 1]);
        for ispeed = 1: length(speeds_sort)
            for coord = 1:3
                subplot(2,3,coord)
                % experimental grf                
                Fsel = Dat(ispeed).benchmark.GRF_r(:,coord); 
                Fsel = Fsel./mexp.*msim;
                plot(Fsel,'Color',Colors_Speeds(ispeed,:));hold on;

                % simulated grf
                subplot(2,3,coord+3)
                dsel = Dat(ispeed).R.ground_reaction.GRF_r(:,coord);
                dsel_int = interp1(1:length(dsel),dsel,linspace(1,length(dsel),100));
                l(ispeed) =plot(dsel_int,'Color',Colors_Speeds(ispeed,:));hold on;


            end
        end
        legend(l,WalkSpeed_legend,'Interpreter','tex');
        legend boxoff;
        title_grf = {'GRFx','GRFy','GRFz'};
        for isubpl =1:6
            subplot(2, 3,isubpl)
            set(gca,'box','off');
            set(gca,'FontSize',10);
            if isubpl<=3
                title(title_grf{isubpl},'interpreter','none');
            else
                xlabel('% gait cycle');
            end
            if isubpl == 1
                ylabel({'experiment','force [N]'});
            elseif isubpl == 4
                ylabel({'simulation','force [N]'});
            end
        end

        % plot stride frequency
        figure('Name','vanderzee: stride frequency','Color',[1 1 1]);

        % experimental stride frequency
        exp_freq = nan(length(speeds),1);
        sim_freq = nan(length(speeds),1);
        for ispeed = 1:length(speeds)
            exp_freq(ispeed) = Dat(ispeed).benchmark.stride_frequency ...
                ./(sqrt(g/Lexp)) .* (sqrt(g/Lsim)); % exp frequency looks weird here compared to other studies
            sim_freq(ispeed) = Dat(ispeed).R.spatiotemp.stride_freq;
        end
        Cs = [0 0 0];
        mk = 4;
        plot([min(exp_freq) max(exp_freq)], [min(exp_freq) max(exp_freq)],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        plot(exp_freq,sim_freq,'ok','Color',Cs,'MarkerFaceColor',Cs,...
            'MarkerSize',mk);        
        set(gca,'box','off');
        set(gca,'FontSize',10);
        xlabel('measured stride frequency');
        ylabel('simulated stride frequency');



        % plot metabolic power

        clear Dat
    end
end




end