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



% The only input is S and osim_path, maybe additional input S_benchmark
S_input = S; % make copy of S
osim_path_input = osim_path;

% save all settings for the benchmark procedure in the output folder. We
% will use this settings once all simulations are finished to benchmark the
% simulation results / compare it to experimental data
if ~isfolder(S_benchmark.out_folder)
    mkdir(S_benchmark.out_folder)
end
save(fullfile(S_benchmark.out_folder,'benchmark_settings.mat'),'S',...
    'S_benchmark','osim_path');

% get path information

% flow of the function
%-----------------------
% 0. pre-processing default model
%       I think the most clean way is to add additional settings to control
%       the flow in the run_pred_sim script
S.pre_processing_only = True; 
run_pred_sim(S,osim_path);
S.pre_processing_only = False;

%
% 1. check if we have to make multiple .dll files
%    (this is needed for walking on a slope and with added mass)
%    if this is the case:
%       - make new .osim files
%       - create new .dll files
%       - copy preprocessing information from default model to new models
%           - F_Falisse_et_al_2022_IO.mat
%           - f_lMT_vMT_dM_poly_3_9

% 2. Run batch simulations

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
        for ispeed = 1:length(gait_speeds)
            % start from default input settings
            S = S_input;
            % set forward velocity
            S.misc.forward_velocity = gait_speeds(i_speed);
            % adapt folder to save results
            S.misc.save_folder  = fullfile(S_benchmark.out_folder,'vanderzee2022',...
                ['speed_' num2str(round(gait_speeds(i_speed)*100)) '_id_' num2str(id_trials)]);
        end
    end
end

% batch simulation gait speeds
if S_benchmark.gait_speeds
    
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
    for i_speed = 1:length(gait_speeds)  
        % start from default input settings
        S = S_input;
        % set forward velocity
        S.misc.forward_velocity = gait_speeds(i_speed);
        % adapt folder to save results
        S.misc.save_folder  = fullfile(S_benchmark.out_folder,'gait_speeds',...
            ['speed_' num2str(round(gait_speeds(i_speed)*100))]); 
    end
end


end