% --------------------------------------------------------------------------
% ankle_exoskeleton
%   Simulate walking with an ankle exoskeleton that provides a
%   plantarflexion torque in function of the progression within a stride.
%   Running this simulation requires the function desired_torque_generator,
%   which can be downloaded from
%   https://www.science.org/doi/suppl/10.1126/science.aal5054/suppl_file/aal5054_zhang_sm_data_s2.zip.
%   Note that using the downloaded function will throw an error related to
%   the variable prev_stride_time. To fix this, set it equal to
%   stride_period.
%
%   References
%   [1] J. Zhang et al., “Human-in-the-loop optimization of exoskeleton 
%   assistance during walking,” Science, vol. 356, pp. 1280–1283, Jun. 2017, 
%   doi: 10.1126/science.aal5054.
% 
%   See also desired_torque_generator
%
% Original author: Lars D'Hondt
% Original date: 16/September/2024
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
addpath(fullfile(pathRepo,'DefaultSettings'))

[S] = initializeSettings('DHondt_et_al_2024_3seg');

%% Settings

% name of the subject
S.subject.name = 'DHondt_et_al_2024_3seg';

% path to folder where you want to store the results of the OCP
S.misc.save_folder  = fullfile(pathRepoFolder,'PredSimResults',S.subject.name); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.solver.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
S.solver.IG_selection_gaitCyclePercent = 100;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);


%% Add ankle exoskeleton

% select orthosis function
exo1.function_name = 'ankleExoZhang2017';

% set path to downloaded function - CHANGE THIS
exo1.dependencies_path = 'c:/path/to/oaal5054_zhang_sm_data_s2';

% set parameters of assistance profile
exo1.peak_torque = 20; % [Nm]
exo1.peak_time = 52.9; % [%]
exo1.rise_time = 26.2; % [%]
exo1.fall_time = 9.8;  % [%]

% add orthosis on right side
exo1.left_right = 'r';
S.orthosis.settings{1} = exo1;

% add the same orthosis on left side
exo1.left_right = 'l';
S.orthosis.settings{2} = exo1;

%% Run predictive simulations

[savename] = runPredSim(S, osim_path);

