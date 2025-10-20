%% Predictive Simulations of Human Gait

% This script starts the predictive simulation of human movement. The
% required inputs are necessary to start the simulations. Optional inputs,
% if left empty, will be taken from getDefaultSettings.m.

clear
close all
clc

% Check BLAS/LAPACK version; add functions from LinearAlgebra subdirectory
% to path in case Intel is *not* used
blas_version = version('-blas')
lapack_version = version('-lapack')
if ~startsWith(lapack_version, 'Intel')
    addpath(fullfile(getenv('PWD'), 'LinearAlgebra'))
end

% path to the repository folder
[pathRepo,~,~] = fileparts(mfilename('fullpath'));
% path to the folder that contains the repository folder
[pathRepoFolder,~,~] = fileparts(pathRepo);

% if the OpenSim module is loaded, make its Java library available
if isenv('EBROOTOPENSIM')
    javaclasspath(fullfile(getenv('EBROOTOPENSIM'), 'sdk', 'Java', 'org-opensim-modeling.jar'));
end

% if the CasADi-MATLAB module is loaded, expose its matlab bindings
if isenv('EBROOTCASADIMINMATLAB')
    addpath(fullfile(getenv('EBROOTCASADIMINMATLAB'), 'matlab'))
end

%% Initialize S
addpath(fullfile(pathRepo,'DefaultSettings'))

[S] = initializeSettings('Falisse_et_al_2022');

%% Settings (custom, defaults are defined in getDefaultSettings.m)

% name of the subject
S.subject.name = 'Falisse_et_al_2022';

% path to folder where you want to store the results of the OCP
S.misc.save_folder  = fullfile(pathRepoFolder,'PredSimResults',S.subject.name); 

% gait cycle
S.misc.gaitmotion_type = 'FullGaitCycle';
%S.misc.gaitmotion_type = 'HalfGaitCycle';

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
S.solver.IG_selection_gaitCyclePercent = 100;
% S.solver.IG_selection = 'quasi-random';

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% Mesh settings
S.solver.N_meshes = 5;

%% Add exoskeleton

% select orthosis function
exo1.function_name = 'hipexo';

% set parameters
exo1.dynamics.xl = -17; % [Nm]
exo1.dynamics.xu = 17; % [Nm]
exo1.dynamics.ul = -15; % [Nm]
exo1.dynamics.uu = 15; % [Nm]

exo1.gain = 1; %EMG gain for gluteus maximus

% add orthosis on left side
exo1.left_right = 'l';
S.orthosis.settings{1} = exo1;

% add the same orthosis on right side
exo1.left_right = 'r';
S.orthosis.settings{2} = exo1;

%% Run predictive simulations

[savename] = runPredSim(S, osim_path);


%% Plot results
% see .\PlotFigures\run_this_file_to_plot_figures.m for more

% If SLURM_JOB_ID is set, this means we are inside a Slurm job context and the
% default way to plot figures will not work. That is probably fixable, but in
% the meantime, skip the plotting in this case.

if (~S.solver.run_as_batch_job) && (~isenv('SLURM_JOB_ID'))

    % set path to reference result
    result_paths{1} = fullfile(pathRepo,'Tests','ReferenceResults',...
        'Falisse_et_al_2022','Falisse_et_al_2022_paper.mat');
    
    % set path to saved result
    result_paths{2} = fullfile(S.misc.save_folder,[savename '.mat']);
    
    % Cell array with legend name for each result
    legend_names = {'Reference result', 'Your first simulation'};
    
    % add path to subfolder with plotting functions
    addpath(fullfile(pathRepo,'PlotFigures'))
    
    figure_settings(1).name = 'all_angles';
    figure_settings(1).dofs = {'all_coords'};
    figure_settings(1).variables = {'Qs'};
    figure_settings(1).savepath = [];
    figure_settings(1).filetype = {};
    
    % call plotting function
    plot_figures(result_paths, legend_names, figure_settings);

end
