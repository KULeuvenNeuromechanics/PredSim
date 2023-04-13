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
pathDefaultSettings = fullfile(pathRepo,'DefaultSettings');
addpath(pathDefaultSettings)

[S] = initializeSettings();
S.misc.main_path = pathRepo;

addpath(fullfile(S.misc.main_path,'VariousFunctions'))

%% Required inputs
% name of the subject
S.subject.name = 'DHondt_2023_3seg';

% path to folder where you want to store the results of the OCP
S.subject.save_folder  = fullfile(pathRepoFolder,'PredSimResults',S.subject.name); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.subject.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
S.subject.IG_selection_gaitCyclePercent = 100;
% S.subject.IG_selection = 'quasi-random';

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% path to folder with program to create dll files from opensim model (will
% be downloaded automatically if it is not there)
S.Cpp2Dll.PathCpp2Dll_Exe = 'C:\GBW_MyPrograms\Osim2Dll_exe';

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
S.bounds.coordinates        = {'pelvis_ty',0.55,1.1, 'pelvis_tilt',-2.9302,nan};

% % S.metabolicE - metabolic energy
S.metabolicE.tanh_b = 100;
% S.metabolicE.model  = '';

% % S.misc - miscellanious
% S.misc.v_max_s             = ;
% S.misc.visualize_bounds    = 1;
% S.misc.gaitmotion_type     = '';
% S.misc.msk_geom_eq         = '';
% S.misc.poly_order.lower    = ;
% S.misc.poly_order.upper    = ;
% S.misc.msk_geom_bounds     = {{'knee_angle_r','knee_angle_l'},-120,10,'pelvis_tilt',-30,30};
% S.misc.msk_geom_n_samples = ;
% S.misc.gaitmotion_type = 'FullGaitCycle';

% % S.post_process
S.post_process.make_plot = 0;
% S.post_process.savename  = 'datetime';
% S.post_process.load_prev_opti_vars = 1;
% S.post_process.rerun   = 1;
% S.post_process.result_filename = '';

% % S.solver
% S.solver.linear_solver  = '';
% S.solver.tol_ipopt      = ;
% S.solver.max_iter       = 5;
% S.solver.parallel_mode  = '';
% S.solver.N_threads      = 6;
% S.solver.N_meshes       = 100;
% S.solver.par_cluster_name = ;
S.solver.CasADi_path    = 'C:\GBW_MyPrograms\casadi_3_5_5';


% % S.subject
% S.subject.mass              = ;
% S.subject.IG_pelvis_y       = ;
S.subject.v_pelvis_x_trgt   = 1.33;
% S.subject.IK_Bounds = ;
% S.subject.muscle_strength   = ;
S.subject.muscle_pass_stiff_shift = {{'soleus','_gas','per_','tib_','_dig_','_hal_','FDB'},0.9}; %,'FDB'
% S.subject.muscle_pass_stiff_scale = ;
S.subject.tendon_stiff_scale      = {{'soleus','_gas'},0.5};
% S.subject.scale_MT_params = {{'soleus_l'},'FMo',0.9,{'soleus_l'},'alphao',1.1};
% S.subject.spasticity        = ;
% S.subject.muscle_coordination = ;
% S.subject.mtp_type          = '2022paper';
% S.subject.set_stiffness_coefficient_selected_dofs = {'mtp_angle',1};
% S.subject.set_damping_coefficient_selected_dofs = {'mtp_angle',2};
% S.subject.set_limit_torque_coefficients_selected_dofs = {{'mtj_angle_l','mtj_angle_r'},[0,0,0,0],[0,0]};

% % S.weights
% S.weights.E         = ;
% S.weights.E_exp     = ;
% S.weights.q_dotdot  = ;
% S.weights.e_arm     = ;
% S.weights.pass_torq = ;
% S.weights.a         = ;
% S.weights.slack_ctrl = ;
% S.weights.pass_torq_includes_damping = ;

% %S.Cpp2Dll: required inputs to convert .osim to .dll
% S.Cpp2Dll.compiler = 'Visual Studio 17 2022';
% S.Cpp2Dll.export3DSegmentOrigins = ;
S.Cpp2Dll.verbose_mode = 0; % 0 for no outputs from cmake
% S.Cpp2Dll.jointsOrder = {'ground_pelvis', 'hip_l', 'hip_r', 'knee_l', 'knee_r',...
%               'ankle_l', 'ankle_r', 'subtalar_l', 'subtalar_r', 'midtarsal_l',...
%               'midtarsal_r', 'tarsometatarsal_l', 'tarsometatarsal_r',...
%               'mtp_l', 'mtp_r', 'back', 'acromial_l', 'acromial_r', 'elbow_l',...
%               'elbow_r', 'radioulnar_l', 'radioulnar_r', 'radius_hand_l',...
%               'radius_hand_r'};
% S.Cpp2Dll.coordinatesOrder = {'pelvis_tilt', 'pelvis_list', 'pelvis_rotation',...
%     'pelvis_tx', 'pelvis_ty', 'pelvis_tz', 'hip_flexion_l', 'hip_adduction_l',...
%     'hip_rotation_l', 'hip_flexion_r', 'hip_adduction_r', 'hip_rotation_r',...
%     'knee_angle_l','knee_angle_r', 'ankle_angle_l', 'ankle_angle_r',...
%     'subtalar_angle_l','subtalar_angle_r',...
%     'mtj_angle_l', 'mtj_angle_r',...
%     'mtp_angle_l', 'mtp_angle_r',...
%     'lumbar_extension', 'lumbar_bending', 'lumbar_rotation', 'arm_flex_l',...
%     'arm_add_l', 'arm_rot_l', 'arm_flex_r', 'arm_add_r', 'arm_rot_r',...
%     'elbow_flex_l', 'elbow_flex_r'};
        
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

