%% Predictive Simulations of Human Gait

% This script starts the predictive simulation of human movement. The
% required inputs are necessary to start the simulations. Optional inputs,
% if left empty, will be taken from getDefaultSettings.m.

clear
close all
clc

% path to the repository folder
pathRepo ='C:\GBW_MyPrograms\sharedPredSim';
% path to the folder that contains the repository folder
[pathRepoFolder,~,~] = fileparts(pathRepo);

%% Initialize S
[S] = initializeSettings('DHondt_et_al_2024_3seg');

addpath(fullfile(pathRepo,'DefaultSettings'))

%% Settings
% name of the subject
S.subject.name = 'CP_AFO_3';

% path to folder where you want to store the results of the OCP
S.misc.save_folder  = fullfile(pathRepoFolder,'PredSimResults',S.subject.name); 

S.misc.main_path = pathRepo ;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '_PredSim.osim']);
% osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);
% 
S.settings.IG = 'ig1';
S.settings.scaling = 'sc1';
S.settings.jointparams ='j1';
S.settings.set_limit_torque_coefficients_selected_dofs = 'limtor1';
% 
% %%%%%%
% % additional
S.settings.muscle_strength = 'IWA_NOmmt';
S.settings.scale_MT_params = 'mt2';
% S.settings.motor_control = 'SynW';
%%%%
% ---------------------------------------------------------------
settings_path = fullfile(pathRepo,'Subjects',S.subject.name,['settings_',S.subject.name,'.m']);

% S.post_process.make_plot = 1;
S.misc.savename  = 'datetime_job';
S.solver.CasADi_path = 'C:\GBW_MyPrograms\casadi_3_5_5';
S.solver.par_cluster_name = '6x4';
S.solver.run_as_batch_job = 1;
% S.solver.batch_job_paths = {};


S = settings_CP_AFO_3(S);
S.solver.adaptIG = true;

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
