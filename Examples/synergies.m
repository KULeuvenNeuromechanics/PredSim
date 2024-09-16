% --------------------------------------------------------------------------
% Implementation of motor control modelled through muscle synergies
% 
%   Here we provide an example of how to define the settings to run a
%   predictive simulation using muscle synergies to model motor control. 
% 
%   In this specific example, we run a simulation of a half gait cycle (symmetric)
%   with a 2D model, imposing four synergies per leg, and tracking synergy
%   weights. The values for muscle weights have been obtained from the figures 
%   from Pitto et al. 2020, and correspond to the generic weights for four 
%   synergies of a group of children with cerebral palsy.
% 
%   Pitto L, van Rossom S, Desloovere K, Molenaers G, Huenaerts C, et al. (2020) 
%   Pre-treatment EMG can be used to model post-treatment muscle coordination 
%   during walking in children with cerebral palsy. PLOS ONE 15(2): e0228851. 
%   https://doi.org/10.1371/journal.pone.0228851
%
%   The full description of the different settings is given in 
%   Documentation/SettingsOverview.md
% 
% Original author: Míriam Febrer-Nafría
% Original date: 16/Sept/2024
% --------------------------------------------------------------------------

clear
close all
clc

[pathExDir,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathExDir);
[pathRepoFolder,~,~] = fileparts(pathRepo);

addpath(fullfile(pathRepo,'DefaultSettings'))
addpath(pathRepo)

%% Initialize S

[S] = initializeSettings('gait1018');

%% Settings

% name of the subject
S.subject.name = 'gait1018';

% path to folder where you want to store the results of the OCP
S.misc.save_folder  = fullfile(pathExDir,'ExampleResults','Synergies'); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
S.solver.IG_selection_gaitCyclePercent = 100;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);


%% Synergy settings
% Use muscle synergies
S.subject.synergies = 1; 

% Impose 4 synergies for right side (and 4 for left side - symmetric)
S.subject.NSyn_r = 4;

% Track synergy weights
S.subject.TrackSynW = 1; 

% Number of tracked synergies (may be different from the number of synergies)
S.subject.TrackSynW_NSyn_r = 4; 

% Synergy weights to be tracked
% The names of muscles should be updated for each model
S.subject.knownSynW_r = {'rect_fem_r', [0 0.28 0 0.75],...
    'vasti_r', [0.1 0.1 0.05 0.8],...
    'bifemsh_r', [0.6 0.4 0.15 0.2],...
    'hamstrings_r', [0 0 1 0],...
    'tib_ant_r', [0.2 1 0.01 0],...
    'gastroc_r', [0.95 0.08 0.04 0.03],...
    'soleus_r', [1 0.05 0 0.05],...
    'glut_max_r', [0.15 0 0.04 1]};

%% Run predictive simulations

[savename] = runPredSim(S, osim_path);

%% Plot results
% see .\PlotFigures\run_this_file_to_plot_figures.m for more

if ~S.solver.run_as_batch_job

    % set path to reference result
    result_paths{1} = fullfile(pathRepo,'Subjects','gait1018',...
        'gait1018_v1.mat');
    
    % set path to saved result
    result_paths{2} = fullfile(S.misc.save_folder,[savename '.mat']);
    
    % Cell array with legend name for each result
    legend_names = {'Reference result',...
        sprintf('%i muscle synergies',S.subject.NSyn_r)};
    
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
