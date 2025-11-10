%% Example script to benchmark Falisse 2022

% paper: Modeling toes contributes to realistic stance knee mechanics in
% three-dimensional predictive simulations of walking 
% (https://doi.org/10.1371/journal.pone.0256311)

% benchmarking paper: [ref]


%% Inputs: Model definition and general settings

%----------     Path information ----------------------
% path to the repository folder
[pathRepo_temp,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathRepo_temp);
% path to the folder that contains the repository folder
[pathRepoFolder,~,~] = fileparts(pathRepo);
addpath(pathRepo);
addpath(fullfile(pathRepo,'Benchmarking'));

%----------     Model settings ----------------------
% add folder with default settings and initialise settings for Falisse 2022
addpath(fullfile(pathRepo,'DefaultSettings'));

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
S.solver.N_meshes       = 50;


%% Information for batch processing simulations
%----------     Solver information ------------------
S.solver.run_as_batch_job = true;
S.solver.N_threads      = 2;
S.solver.par_cluster_name = 'Cores4';

%% Specific settings for benchmark function

% all these settings will be treated as optional. If not provided I will
% assume that the user wants to benchmark everything (gait speed, slope,
% walking with added mass).

% benchmark specific studies
S_benchmark.studies = {'vanderzee2022','browning2008','koelewijn2019',...
    'gomenuka2014','schertzer2014'};
% options are:
%   vanderzee2022: variations in gait speed
%   koelewijn2019: variation in gait speed and slope
%   browning2008: added mass to body segments
%   schertzer2014: added mass to body segments and variations in gait speed
%   gomenuka2014: added mass to pelvis, walking on a slope and various
%   speeds


% % benchmark gait speed simulations
S_benchmark.gait_speeds = true;
S_benchmark.gait_speed_range = [0.6 2];
S_benchmark.gait_speeds_selection = 0.6:0.2:2;
% 
% % benchmark walking on a slope
% S_benchmark.slope = false;
% S_benchmark.slope_range = [0 24]; % in %
% S_benchmark.slope_range_selection = 0:2:24;
% 
% % benchmarking walking with added mass
% S_benchmark.added_mass = false;
% S_benchmark.added_mass_studies = {'Browning','Schertzer','Huang'}; % add studies you want to replicate here

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


%% Run benchmarking procedure

% matlab function used to start all simulations
benchmark_predsim(S,osim_path,S_benchmark);


%% Compare simulation with experiment

% you can use the function benchmark_results to do this. This function
% adds the experimental data to your simulations and has some plot
% utilities to compare simulations and experiments. 
% please run this function when all simulations are finished

%benchmark_results(S_benchmark.out_folder,'BoolPlot',true);