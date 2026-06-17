%% Predictive Simulations of Human Gait

% This script starts the predictive simulation of human movement. The
% required inputs are necessary to start the simulations. Optional inputs,
% if left empty, will be taken from getDefaultSettings.m.

clear
close all
clc
tic
% BAK CHANGES 
% copied stack overflow
cmdir = 'C:\Program Files\CMake\bin\';
setenv('PATH', [cmdir ';' getenv('PATH')]);
% OS - I HAVE PUT THIS IN MY SYS ENV
osimdir = 'C:\OpenSim 4.4\bin\';
setenv('PATH', [osimdir ';' getenv('PATH')]);


    % path to the repository folder
[pathRepo,~,~] = fileparts(mfilename('fullpath'));
% path to the folder that contains the repository folder
[pathRepoFolder,~,~] = fileparts(pathRepo);

%% Initialize S
pathDefaultSettings = fullfile(pathRepo,'DefaultSettings');
addpath(pathDefaultSettings)

[S] = initializeSettings();
S.misc.main_path = pathRepo;

addpath(fullfile(S.misc.main_path,'VariousFunctions'))

%% Required inputs
% name of the subject
S.subject.name = 'Falisse_et_al_2022_gaitRetraining';

% path to folder where you want to store the results of the OCP
S.subject.save_folder  = fullfile(pathRepoFolder,'PredSimResults',S.subject.name); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess

% Typiclaly used a predsim solution using
S.subject.IG_selection = "quasi-random";

S.subject.IG_selection_gaitCyclePercent = 100;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% path to folder with program to create dll files from opensim model (will
% be downloaded automatically if it is not there)
S.Cpp2Dll.PathCpp2Dll_Exe = 'C:\Users\YourUserName\Documents\MATLAB\osimModelToExtFunc';

% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = 1;

%% Optional inputs
% see README.md in the main folder for information about these optional
% inputs.

% % S.bounds
% S.bounds.a.lower            = ;
% S.bounds.calcn_dist.lower   = ;
% S.bounds.toes_dist.lower    = ;
% S.bounds.tibia_dist.lower   = ;
% S.bounds.SLL.upper          = ;
% S.bounds.SLR.upper          = ;
% S.bounds.dist_trav.lower    = ;
% S.bounds.t_final.upper      = ;
% S.bounds.t_final.lower      = ;
% S.bounds.coordinates        = {'pelvis_tilt',-30,30,'pelvis_list',-30,30};

%Added for gait retraining implementation
S.bounds.coordinates        = { 'medial_knee_ty_l',-0.001,0.001,  'medial_knee_ty_r' , -0.001,0.001, ...
                               'int_knee_ty_l',-0.001,0.001,  'int_knee_ty_r' , -0.001,0.001, ...
                               'subtalar_angle_r',-0.001,0.001,  'subtalar_angle_l' , -0.001,0.001}; 


% % S.metabolicE - metabolic energy
% S.metabolicE.tanh_b = ;
% S.metabolicE.model  = '';

% % S.misc - miscellanious
% S.misc.v_max_s             = ;
% S.misc.visualize_bounds    = 1;
S.misc.gaitmotion_type     = 'FullGaitCycle'; % Half or Full
% S.misc.msk_geom_eq         = '';
% S.misc.poly_order.lower    = 1;
S.misc.poly_order.upper    = 10;
S.misc.msk_geom_bounds      = {{'knee_angle_r','knee_angle_l'},-120,10,{'hip_flexion_r', 'hip_flexion_l'} ,-90,90 ...
    {'ankle_angle_r' , 'ankle_angle_l'}, -30,30, {'int_knee_ty_l', 'int_knee_ty_r'}, -0.01, 0.01 ...
    {'medial_knee_ty_l', 'medial_knee_ty_r'}, -0.01, 0.01};

% %Added for gait retraining implementation
S.misc.jcf_to_min             = {'medial_knee_ty_r','medial_knee_ty_l','int_knee_ty_r','int_knee_ty_l'};
S.misc.jcf_to_min_weights     = struct('medial_knee_ty_r',1,'medial_knee_ty_l',1, 'int_knee_ty_r',1, 'int_knee_ty_l',1) ;

% adding toggle for tracking similar to Lars imp.
S.TrackSim                      =  1; 
% set to as IG for now.

S.Track.Qref                    = 'Directory/to/your/subjects/initial/habitual/gait/pattern.mot';

                                                                                                          
% %Added for gait retraining implementation
%S.Track.dofs                    = {'hip_flexion_r',  'ankle_angle_r','knee_angle_r','hip_flexion_l',  'ankle_angle_l','knee_angle_l' };
S.Track.dofs                    = {'hip_flexion_r'};

% Optianl weights - if not set = 20 
%S.Track.dofs_weights            = struct('hip_flexion_r',5,'ankle_angle_r',20,'knee_angle_r',5,'hip_flexion_l',5, 'ankle_angle_l',20,'knee_angle_l', 5);
S.Track.dofs_weights            = struct('hip_flexion_r',0);
%%
% this takes alot fo space to save tbh 
%S.misc.saveWS_for_rerun = 0;

% % S.post_process
S.post_process.make_plot = 0;
S.post_process.savename  = 'datetime';

% use if doesnt output kinematics 
%S.post_process.rerun   = 1;
%S.post_process.result_filename = 'Falisse_et_al_2022_20240822T161534';

% % S.solver
% S.solver.linear_solver  = '';
% S.solver.tol_ipopt      = ;
S.solver.max_iter       = 500 ;
% S.solver.parallel_mode  = ''; 
% S.solver.N_threads      = 6;
% S.solver.N_meshes       = ;
% S.solver.par_cluster_name = ;
S.solver.CasADi_path    = 'C:\Users\YourUserName\Documents\MATLAB\Casadi';


% % S.subject
 S.subject.mass              = 75;
% S.subject.IG_pelvis_y       = ;
S.subject.v_pelvis_x_trgt   = 1.12; % CHANGED BASED ON INPUT IK 1.33;

% S.subject.IK_Bounds = ;
% S.subject.muscle_strength   = ;
% S.subject.muscle_pass_stiff_shift = {{'soleus_l','soleus_r'},0.9,{'tib_ant_l'},1.1};
% S.subject.muscle_pass_stiff_scale = ;
% S.subject.tendon_stiff_scale      = ;
S.subject.mtp_type          = '2022paper';
% S.subject.scale_MT_params         = {{'soleus_l'},'FMo',0.9,{'soleus_l'},'alphao',1.1};
% S.subject.spasticity        = ;
% S.subject.muscle_coordination = ;
S.subject.set_stiffness_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},25};
S.subject.set_damping_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},2};
% S.subject.set_limit_torque_coefficients_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},[0,0,0,0],[0,0]};

% % S.weight s           ; THESE ARE THE DEAULTS VALUES 
%S.weights.E         = 250; %500
%S.weights.E_exp     = 1;  %2
%S.weights.q_dotdot  = 25000; %50000
%S.weights.e_arm     = 0; %10000000
%S.weights.pass_torq = 0; %1000
%S.weights.a         = 20; %2000
%S.weights.slack_ctrl = 0 ; %0.0001
% S.weights.pass_torq_includes_damping = ;

% %S.Cpp2Dll: required inputs to convert .osim to .dll

%S.Cpp2Dll.compiler = 'Visual Studio 17 2022';
%S.Cpp2Dll.compiler = 'Microsoft Visual Studio 12.0';

% S.Cpp2Dll.export3DSegmentOrigins = ;
S.Cpp2Dll.verbose_mode = 0; % 0 for no outputs from cmake
% S.Cpp2Dll.jointsOrder = ;
% S.Cpp2Dll.coordinatesOrder = ;
        
%% Run predictive simulations
% Check for updates in osim2dll
S.Cpp2Dll.PathCpp2Dll_Exe = InstallOsim2Dll_Exe(S.Cpp2Dll.PathCpp2Dll_Exe);


% Start simulation
if S.solver.run_as_batch_job
    add_pred_sim_to_batch(S,osim_path)
else
    [savename] = run_pred_sim(S,osim_path);
end

%% Plot results
if S.post_process.make_plot && ~S.solver.run_as_batch_job
    % set path to saved result
    result_paths{2} = fullfile(S.subject.save_folder,[savename '.mat']);
    % add path to subfolder with plotting functions
    addpath(fullfile(S.misc.main_path,'PlotFigures'))
    % call plotting script
    run_this_file_to_plot_figures
end

toc