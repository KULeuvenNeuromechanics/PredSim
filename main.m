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
pathDefaultSettings = [pathRepo '\DefaultSettings'];
addpath(pathDefaultSettings)

[S] = initializeSettings();
S.misc.main_path = pathRepo;

addpath([S.misc.main_path '\VariousFunctions'])

%% Required inputs
% name of the subject
S.subject.name = 'Subject1_2D';

% path to folder where you want to store the results of the OCP
S.subject.save_folder  = fullfile(pathRepoFolder,'PredSimResults',S.subject.name); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.subject.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Default.mot');
S.subject.IG_selection_gaitCyclePercent = 50;


% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = 0;

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

% % S.metabolicE - metabolic energy
% S.metabolicE.tanh_b = ;
% S.metabolicE.model  = '';

% % S.misc - miscellanious
% S.misc.v_max_s             = ;
% S.misc.visualize_bounds    = 1;
% S.misc.gaitmotion_type     = 'FullGaitCycle';
% S.misc.msk_geom_eq         = '';
% S.misc.poly_order.lower    = ;
% S.misc.poly_order.upper    = ;
% S.misc.msk_geom_bounds      = {{'knee_angle_r','knee_angle_l'},-120,10,'pelvis_tilt',-30,30};
% S.misc.coordPosZeroForExternal = {'pelvis_tx'};

% % S.post_process
% S.post_process.make_plot = '';
% S.post_process.savename  = 'datetime';
% S.post_process.rerun   = 1;
% S.post_process.rerun_from_w = 1;
S.post_process.result_filename = 'Subject1_2D_jit_Ox_p6_openmp';
% S.post_process.result_filename = 'Subject1_2D_no_jit_p6';

% % S.solver
% S.solver.linear_solver  = 'ma86';
% S.solver.tol_ipopt      = ;
% S.solver.max_iter       = 5;
% S.solver.parallel_mode  = '';
S.solver.N_threads      = 6;
% S.solver.N_meshes       = 50;
% S.solver.par_cluster_name = ;
S.solver.CasADi_path    = 'C:\GBW_MyPrograms\casadi_3_5_5';
S.Solver.jit = 1;

% % S.subject
% S.subject.mass              = ;
% S.subject.IG_pelvis_y       = ;
S.subject.v_pelvis_x_trgt   = 1.33;
% S.subject.IK_Bounds = ;
% S.subject.muscle_strength   = ;
% S.subject.muscle_pass_stiff_shift = {{'soleus_l','soleus_r'},0.9,{'tib_ant_l'},1.1};
% S.subject.muscle_pass_stiff_scale = ;
% S.subject.tendon_stiff_scale      = {{'soleus_l','soleus_r','gastroc_r','gastroc_l'},0.5};
% S.subject.mtp_type          = '2022paper';
% S.subject.MT_params         = ;
% S.subject.spasticity        = ;
% S.subject.muscle_coordination = ;
% S.subject.set_stiffness_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},25};
% S.subject.set_damping_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},2};
% S.subject.set_limit_torque_coefficients_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},[0,0,0,0],[0,0]};

% % S.weights
% S.weights.E         = 0;
% S.weights.E_exp     = 1;
% S.weights.q_dotdot  = 1;
% S.weights.e_arm     = 1;
% S.weights.pass_torq = 0;
% S.weights.a         = 1;
% S.weights.slack_ctrl = ;

% %S.Cpp2Dll: required inputs to convert .osim to .dll
% optional: if you want to install the opensimExe
S.Cpp2Dll.PathCpp2Dll_Exe = InstallOsim2Dll_Exe('C:\GBW_MyPrograms\Osim2Dll_exe'); %(optional: if you want to install the opensimExe)
% S.Cpp2Dll.compiler = 'Visual Studio 15 2017 Win64';
S.Cpp2Dll.export3DSegmentOrigins = [];
S.Cpp2Dll.verbose_mode = 0; % 0 for no outputs from cmake

%% Run predictive simulations
if S.solver.run_as_batch_job
    add_pred_sim_to_batch(S,osim_path)
else
    run_pred_sim(S,osim_path);
end


