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
S_benchmark.studies = {'vanderzee2022'};
% options are:
%   vanderzee2022: variations in gait speed
%   koelewijn2019: variation in gait speed and slope
%   browning: added mass to body segments
%   schertzer: added mass to body segments and variations in gait speed


% benchmark gait speed simulations
S_benchmark.gait_speeds = false;
S_benchmark.gait_speed_range = [0.6 2];
S_benchmark.gait_speeds_selection = 0.6:0.2:2;

% benchmark walking on a slope
S_benchmark.slope = false;
S_benchmark.slope_range = [0 24]; % in %
S_benchmark.slope_range_selection = 0:2:24;

% benchmarking walking with added mass
S_benchmark.added_mass = false;
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


%% Run benchmarking procedure

% matlab function used to start all simulations
benchmark_predsim(S,osim_path,S_benchmark);


%% Compare simulation with experiment




% 3. method to attach data from literature to each simulation or make a
% table as in my script TableSim.m (benchmarking study). Not sure how to do
% this properply.
%
% Other option is to only provide a graphical benchmarking. In other words
% only make a figure with measured and simulated values. For each
% experimental condition look for the datapoint that is most similar to 
% measured values
%
% current idea is to write function here that attaches experimental data to
% simulation results structure. This is only for the replication of
% specific studies. I you want to benchmark generic effects of gait speed I
% will only add the figure (e.g. measured and simulated COT as a function
% of gait speed).
% this is a seperate function that you should run when all simulations are
% finished

% step 0. download data from Benchmarking paper 2025:
% DOI 10.5281/zenodo.14524619
data_exp_path = fullfile(S_benchmark.out_folder,'exp_data');
if ~isfolder(data_exp_path)
    mkdir(data_exp_path);
end

% general data
LinkZipFile = 'https://zenodo.org/records/14524620/files/PredSim_gait_conditions-master_19_12_2024.zip?download=1';
websave(fullfile(data_exp_path,'githubrepo_data.zip'),LinkZipFile);
unzip(fullfile(data_exp_path,'githubrepo_data.zip'),fullfile(data_exp_path,'github_repo_temp'));
delete(fullfile(data_exp_path,'githubrepo_data.zip'));
%       change data format a bit here, we don't need everything
copyfile(fullfile(data_exp_path,'github_repo_temp','PredSim_gait_conditions-master','ExperimentalData'),...
    fullfile(data_exp_path,'general_exp_data'))
rmdir(fullfile(data_exp_path,'github_repo_temp'),'s');

% Van Der Zee kinematic data
if isfield(S_benchmark,'studies') && ~isempty(S_benchmark.studies)
    if any(strcmp(S_benchmark.studies,'vanderzee2022'))
        LinkZipFile = 'https://zenodo.org/records/14524620/files/VanDerZee_Averages_Addbiomech.zip?download=1';
        websave(fullfile(data_exp_path,'vanderzee2022.zip'),LinkZipFile);
        unzip(fullfile(data_exp_path,'vanderzee2022.zip'),fullfile(data_exp_path,'vanderzee2022'));
        delete(fullfile(data_exp_path,'vanderzee2022.zip'));
    end
    % attach experimental data to simulation results
    
end







% 4. (optional) make figures as in 2022 paper (or some of the figures)




