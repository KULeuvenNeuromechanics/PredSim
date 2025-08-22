%% Example script to benchmark Falisse 2022

% paper: Modeling toes contributes to realistic stance knee mechanics in
% three-dimensional predictive simulations of walking 
% (https://doi.org/10.1371/journal.pone.0256311)

% benchmarking paper: [ref]


%% Inputs: Model definition and general settings

%----------     Path information ----------------------
% path to the repository folder
[pathRepo,~,~] = fileparts(mfilename('fullpath'));
% path to the folder that contains the repository folder
[pathRepoFolder,~,~] = fileparts(pathRepo);

%----------     Model settings ----------------------
% add folder with default settings and initialise settings for Falisse 2022
addpath(fullfile(pathRepo,'DefaultSettings'))
[S] = initializeSettings('Falisse_et_al_2022');
% model name
S.subject.name = 'Falisse_et_al_2022';
% path to opensim model
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);
% adapt lower bound on muscle activation (Afschrift 2025)
S.bounds.activation_all_muscles.lower = 0.01;

%----------     Initial guess settings --------------
S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
S.solver.IG_selection_gaitCyclePercent = 100;

%----------     Collocation -------------------------
S.solver.N_meshes       = 100;


%% Information for batch processing simulations
%----------     Solver information ------------------
S.solver.run_as_batch_job = True;
S.solver.N_threads      = 2;
S.solver.par_cluster_name = 'Cores4';

%% Specific settings for benchmark function

% all these settings will be treated as optional. If not provided I will
% assume that the user wants to benchmark everything (gait speed, slope,
% walking with added mass).

% benchmark gait speed simulations
S_benchmark.gait_speeds = True;
S_benchmark.gait_speed_range = [0.6 2];
S_benchmark.gait_speeds_selection = 0.6:0.2:2;

% benchmark walking on a slope
S_benchmark.slope = True;
S_benchmark.slope_range = [0 24]; % in %
S_benchmark.slope_range_selection = 0:2:24;

% benchmarking walking with added mass
S_benchmark.added_mass = True;
S_benchmark.added_mass_studies = {'Browning','Schertzer','Huang'}; % add studies you want to replicate here

% option to replicate specific studies ?. THis is a bit harder as I have to
% compile specific models for each study (e.g. schertzer has walking at
% specific slope with added mass to specific body part)
%   - Koelewijn
%   - Schertzer
%   - Browning
%   - Abe2015
%   - Van Der Zee (gait speeds)


% path information
S_benchmark.out_folder = fullfile(pathRepo,'Results','Benchmark_Falisse2022');


%% Here I will start the function benchmarking

% The only input is S, maybe additional input S_benchmark
S_input = S; % make copy of S

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

    end




end


% 3. method to attach data from literature to each simulation or make a
% table as in my script TableSim.m (benchmarking study). Not sure how to do
% this properply.
%
% Other option to only provide a graphical benchmarking. In other words
% only make a figure with measured and simulated values. For each
% experimental condition look for the datapoint that is most similar to 
% measured values

% 4. (optional) make figures as in 2022 paper (or some of the figures)




