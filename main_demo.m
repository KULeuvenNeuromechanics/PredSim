%% Predictive Simulations of Human Gait

% This script starts the predictive simulation of human movement. The 
% required inputs are necessary to start the simulations. Optional 
% inputs,if left empty, will be taken from getDefaultSettings.m.

clear
close all
clc
% path to the repository folder
% [pathRepo,~,~] = fileparts(mfilename('fullpath'));
pathRepo = 'C:\GBW_MyPrograms\PredSim';
% path to the folder that contains the repository folder
% [pathRepoFolder,~,~] = fileparts(pathRepo);
pathRepoFolder = 'C:\GBW_MyPrograms';

%% Initialize S
pathDefaultSettings = fullfile(pathRepo,'DefaultSettings');
addpath(pathDefaultSettings)

[S] = initializeSettings();
S.misc.main_path = pathRepo;

addpath(fullfile(S.misc.main_path,'VariousFunctions'))

%% Required inputs
% name of the subject
S.subject.name = 'Vitruvian_Man';

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,...
    [S.subject.name '.osim']);

% either choose "quasi-random" or give the path to a .mot file 
% you want to use as initial guess
S.subject.IG_selection = fullfile(S.misc.main_path,'OCP',...
    'IK_Guess_Full_GC.mot');
S.subject.IG_selection_gaitCyclePercent = 100; % 1 gait cycle

% path to folder where you want to store the results of the OCP
S.subject.save_folder = fullfile(pathRepo,'Subjects',S.subject.name); 

% Do you want to run the simulation as a batch job 
% (parallel computing toolbox)
S.solver.run_as_batch_job = 0;

%% Optional inputs
% See README.md in the main folder for information about these 
% optional inputs.

% % Task constraints
% Target forward velocity
S.subject.v_pelvis_x_trgt   = 1.2;

% Cyclic and symmetric
S.misc.gaitmotion_type = 'HalfGaitCycle';


% % Objective function
% Muscle activation exponent
S.weights.a_exp = 2;

% Muscle activation weight
S.weights.a = 180;

% Torque actuator excitation weight
S.weights.e_torqAct = 10;


% % Other
% Filename of saved results
S.post_process.result_filename = 'Demo_TGCS';

% Path to CasADi installation
S.solver.CasADi_path = 'C:\GBW_MyPrograms\casadi_3_5_5';

% Solver can use 2 CPU threads
S.solver.N_threads = 2;








% % Musculoskeletal model
% Set lower bound on muscle activation
S.bounds.activation_all_muscles.lower = 0.01;

% Set coordinate bounds
S.bounds.Qs = {'lumbar_extension',-10,20}; % [°]

% No viscous damping in joints
S.subject.damping_coefficient_all_dofs = 0; 

% No coordinate limit torques
S.subject.set_limit_torque_coefficients_selected_dofs =...
    {'all',[0,0,0,0],[0,0]};

% Do not use these objective function terms
S.weights.E = 0;
S.weights.q_dotdot = 0;
S.weights.pass_torq = 0;

% Do not print outputs from OpenSimAD
S.OpenSimADOptions.verbose_mode = false;


%% Run predictive simulations

% Start simulation
if S.solver.run_as_batch_job
    add_pred_sim_to_batch(S,osim_path)
else
    [savename] = run_pred_sim(S,osim_path);
end


