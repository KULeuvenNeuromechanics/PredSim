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

[S] = initializeSettings('DHondt_et_al_2025');

%% Settings

% name of the subject
S.subject.name = 'DHondt_et_al_2025';

% path to folder where you want to store the results of the OCP
S.misc.save_folder  = fullfile(pathRepoFolder,'PredSimResults',S.subject.name); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
S.solver.IG_selection_gaitCyclePercent = 100;
% S.solver.IG_selection = 'quasi-random';

% S.solver.IG_selection = fullfile(pathRepoFolder,'PredSimResults',...
%     'DHondt_et_al_2024_4seg','DHondt_et_al_2024_4seg_v16.mot');
% S.solver.IG_selection_gaitCyclePercent = 200;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,...
    [S.subject.name '_v4_alt1.osim']);


S.OpenSimADOptions.useSerialisedFunction = true;
S.solver.nlpsol_options.expand = true;

% S.solver.N_threads = 10;

S.solver.ipopt_options.hsllib = ['C:/GBW_MyPrograms/coin-or/coinhsl/build/' ...
        'with-openmp-v2/bin/libhsl.dll'];
S.solver.linear_solver = 'ma97';

% S.weights.E = 0;
% S.solver.N_meshes = 2;
% S.misc.msk_geom_n_samples = 100;
% S.misc.poly_order.upper = 7;

% S.subject.default_coord_lim_torq_coeff = '';

%% Run predictive simulations

[savename] = runPredSim(S, osim_path);


%% Plot results
% see .\PlotFigures\run_this_file_to_plot_figures.m for more

% If SLURM_JOB_ID is set, this means we are inside a Slurm job context and the
% default way to plot figures will not work. That is probably fixable, but in
% the meantime, skip the plotting in this case.

if (~S.solver.run_as_batch_job)

    % set path to reference result
%     result_paths{1} = fullfile(pathRepo,'Tests','ReferenceResults',...
%         'Falisse_et_al_2022','Falisse_et_al_2022_paper.mat');
    
    % set path to saved result
    result_paths{1} = fullfile(S.misc.save_folder,[savename '.mat']);
    
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