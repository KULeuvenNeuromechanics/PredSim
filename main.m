%% Predictive Simulations of Human Gait

% This script starts the predictive simulation of human movement. The
% required inputs are necessary to start the simulations. Optional inputs,
% if left empty, will be taken from getDefaultSettings.m.

clear
close all
clc
% path to the repository folder
[pathRepo,~,~] = fileparts(mfilename('fullpath'));
% path to the folder that contains the repository folder
[pathRepoFolder,~,~] = fileparts(pathRepo);

%% Initialize S
addpath(fullfile(pathRepo,'DefaultSettings'))

[S] = initializeSettings('Falisse_et_al_2022');

%% Settings

% name of the subject
S.subject.name = 'Falisse_et_al_2022';

% path to folder where you want to store the results of the OCP
S.misc.save_folder  = fullfile(pathRepoFolder,'PredSimResults',S.subject.name); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
S.solver.IG_selection_gaitCyclePercent = 100;
% S.solver.IG_selection = 'quasi-random';

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);


%% Run predictive simulations

[savename] = runPredSim(S, osim_path);


%% Plot results
% see .\PlotFigures\run_this_file_to_plot_figures.m for more

if ~S.solver.run_as_batch_job

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
