% --------------------------------------------------------------------------
% synergies
%   Explanation of what this function does. Length depends on function 
%   complexity. Include relevant citations. If applicable, provide an
%   example and/or refer to a unit test. 
%
% 
% Original author: (First name Last name)
% Original date: (Using "30/May/2022" format avoids confusion)
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

[S] = initializeSettings('Falisse_et_al_2022');

%% Settings

% name of the subject
S.subject.name = 'Falisse_et_al_2022';

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
% (from the case of 4 synergies in Pitto et al 2020)
S.subject.knownSynW_r = {'rect_fem_r', [0 0.28 0 0.75],...
    'vas_lat_r', [0.1 0.1 0.05 0.8],...
    {'bifemlh_r','bifemsh_r'}, [0.6 0.4 0.15 0.2],...
    {'semiten_r','semimem_r'}, [0 0 1 0],...
    'tib_ant_r', [0.2 1 0.01 0],...
    {'med_gas_r','lat_gas_r'}, [0.95 0.08 0.04 0.03],...
    'soleus_r', [1 0.05 0 0.05],...
    {'glut_med1_r','glut_med2_r','glut_med3_r'}, [0.15 0 0.04 1]};




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
