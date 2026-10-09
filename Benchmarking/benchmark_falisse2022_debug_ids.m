%% Debug experimental/simulation ID matching for Falisse 2022

% paper: Modeling toes contributes to realistic stance knee mechanics in
% three-dimensional predictive simulations of walking 
% (https://doi.org/10.1371/journal.pone.0256311)

% Follow the same settings as benchmark_falisse2022.m. This script only
% prepares conditions and audits IDs; it never runs or preprocesses simulations.
% Requires MATLAB, but no OpenSim, CasADi or Parallel Computing Toolbox.
% Inspect id_matches in the Variable Editor for side-by-side conditions.
% experimental_coverage lists unused, ambiguous and reused experimental data.
% See Benchmarking/README_ID_DEBUG.md for interpretation and load conventions.


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

% Initialize only what the audit needs; initializeSettings also fetches Git
% refs, which is unnecessary for a local ID check.
S = struct;
S.misc.main_path = pathRepo;
run(fullfile(pathRepo,'Subjects','Falisse_et_al_2022','settings_Falisse_et_al_2022.m'));
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
S.solver.par_cluster_name = 'Cores3'; % use the default local MATLAB cluster

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

% path information
S_benchmark.out_folder = fullfile(pathRepo,'Results','Benchmark_Falisse2022');

% set verbose mode to true
S.OpenSimADOptions.verbose_mode = true;

%% Debug data and report settings
% Leave empty to download/use the standard benchmark data cache. Set this to
% a local JSON dataset folder to work offline or audit a different data copy.
data_folder = '';
% Reports use a separate folder; existing simulation results are never edited.
report_folder = fullfile(pathRepo,'Results','Benchmark_Falisse2022_ID_Debug');

%% Trace benchmarking assignments and compare conditions

[id_matches,experimental_coverage,planned_simulations] = debug_benchmark_ids( ...
    S,osim_path,S_benchmark,'DataFolder',data_folder,'ReportFolder',report_folder);

% Useful filters in the MATLAB Command Window:
% id_matches(id_matches.Status == "condition_mismatch",:)
% id_matches(id_matches.Status == "missing" | id_matches.Status == "ambiguous",:)
% experimental_coverage(experimental_coverage.SelectedStudy & ...
%     experimental_coverage.UniqueMatchCount == 0,:)


%% export analysis results
% export id_matches table to excel

if ~exist(report_folder,'dir')
    mkdir(report_folder);
end

analysis_file = fullfile(report_folder,'falisse2022_id_debug_results.xlsx');

if isfile(analysis_file)
    delete(analysis_file);
end

writetable(id_matches,analysis_file,'Sheet','id_matches');
writetable(experimental_coverage,analysis_file,'Sheet','experimental_coverage');

if isstruct(planned_simulations)
    planned_simulations = struct2table(planned_simulations);
end

writetable(planned_simulations,analysis_file,'Sheet','planned_simulations');