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
S.subject.name = 'Vitruvian_Man_v1';

% path to folder where you want to store the results of the OCP
S.subject.save_folder  = fullfile(pathRepoFolder,'PredSimResults',S.subject.name); 
S.subject.save_folder  = fullfile(S.misc.main_path,'Subjects',S.subject.name,'IG'); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.subject.IG_selection = 'quasi-random';
% S.subject.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
% S.subject.IG_selection_gaitCyclePercent = 100;
S.subject.IG_selection = fullfile(S.misc.main_path,'Subjects','Vitruvian_Man_v1','IG','IG_v25ms_ATx70.mot');
S.subject.IG_selection_gaitCyclePercent = 200;
% S.subject.IG_selection = fullfile('C:\GBW_MyPrograms\OpenSim 4.3\Resources\Code\Matlab\Moco\example2DWalking','referenceCoordinates.mot');
% S.subject.IG_selection_gaitCyclePercent = 50;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);
% S.bounds.activation_all_muscles.lower = 0.01;
% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = 0;

%% Optional inputs
% see README.md in the main folder for information about these optional
% inputs.

% % S.bounds
S.bounds.activation_all_muscles.lower = 0.01;
% S.bounds.SLL.upper          = ;
% S.bounds.SLR.upper          = ;
% S.bounds.dist_trav.lower    = ;
% S.bounds.t_final.upper      = ;
% S.bounds.t_final.lower      = ;
% S.bounds.Qs = {'lumbar_extension',-10,20};

S.subject.v_pelvis_x_trgt   = 3;

S.bounds.Qs = {
%     'pelvis_tilt',-50,20,...
    'lumbar_extension',-20,5,... % v=2
    {'arm_flex_r','arm_flex_l'},-60,60
    };

S.bounds.Qdots = {
    'pelvis_tx',0.01,max(S.subject.v_pelvis_x_trgt)*2,...
%     'pelvis_tilt',-100,100 % v=4.5
    };
    
% S.bounds.Qdotdots = {'pelvis_tilt',-2000,2000};%v=4.5;

S.misc.scaling_Moments = {'all',1.2};

% % to prevent body segments from clipping into eachother
% S.bounds.distanceConstraints(1).point1 = 'calcn_r';
% S.bounds.distanceConstraints(1).point2 = 'calcn_l';
% S.bounds.distanceConstraints(1).direction = 'xz';
% S.bounds.distanceConstraints(1).lower_bound = 0.09;
% S.bounds.distanceConstraints(1).upper_bound = 2;
% 
% S.bounds.distanceConstraints(2).point1 = 'hand_r';
% S.bounds.distanceConstraints(2).point2 = 'femur_r';
% S.bounds.distanceConstraints(2).direction = 'xz';
% S.bounds.distanceConstraints(2).lower_bound = 0.18;
% S.bounds.distanceConstraints(2).upper_bound = 2;
% 
% S.bounds.distanceConstraints(3).point1 = 'hand_l';
% S.bounds.distanceConstraints(3).point2 = 'femur_l';
% S.bounds.distanceConstraints(3).direction = 'xz';
% S.bounds.distanceConstraints(3).lower_bound = 0.18;
% S.bounds.distanceConstraints(3).upper_bound = 2;
% 
% S.bounds.distanceConstraints(4).point1 = 'tibia_r';
% S.bounds.distanceConstraints(4).point2 = 'tibia_l';
% S.bounds.distanceConstraints(4).direction = 'xz';
% S.bounds.distanceConstraints(4).lower_bound = 0.11;
% S.bounds.distanceConstraints(4).upper_bound = 2;
% 
% S.bounds.distanceConstraints(5).point1 = 'toes_r';
% S.bounds.distanceConstraints(5).point2 = 'toes_l';
% S.bounds.distanceConstraints(5).direction = 'xz';
% S.bounds.distanceConstraints(5).lower_bound = 0.1;
% S.bounds.distanceConstraints(5).upper_bound = 2;


% % S.metabolicE - metabolic energy
% S.metabolicE.tanh_b = ;
% S.metabolicE.model  = '';

% % S.misc - miscellanious
% S.misc.v_max_s             = ;
% S.misc.visualize_bounds    = 1;
% S.misc.gaitmotion_type     = '';
% S.misc.msk_geom_eq         = '';
% S.misc.poly_order.lower    = ;
% S.misc.poly_order.upper    = ;
% S.misc.msk_geom_bounds      = {{'knee_angle_r'},0,90,{'mtp_angle_'},-50,20};
% S.misc.default_msk_geom_bounds = fullfile(pathRepo,'Subjects',S.subject.name,'msk_geom_bounds.csv');
% S.misc.gaitmotion_type = 'FullGaitCycle';

% % S.post_process
S.post_process.make_plot = 0;
% S.post_process.savename  = 'datetime';
% S.post_process.load_prev_opti_vars = 1;
% S.post_process.rerun   = 1;
S.post_process.result_filename = ['IG_v' num2str(S.subject.v_pelvis_x_trgt*10) 'ms_ATx70'];

% % S.solver
% S.solver.linear_solver  = '';
S.solver.tol_ipopt      = 3;
% S.solver.max_iter       = 5;
% S.solver.parallel_mode  = '';
% S.solver.N_threads      = 2;
% S.solver.N_meshes       = 100;
% S.solver.par_cluster_name = ;
% S.solver.CasADi_path    = 'C:\GBW_MyPrograms\casadi_3_5_5';


% % S.subject
% S.subject.mass              = ;
% S.subject.IG_pelvis_y       = ;
S.subject.adapt_IG_pelvis_y = 1;
% S.subject.v_pelvis_x_trgt   = 3;
% S.subject.IK_Bounds = ;
% S.subject.muscle_strength   = ;
% S.subject.muscle_pass_stiff_shift = {{'soleus_l','soleus_r'},0.9,{'tib_ant_l'},1.1};
% S.subject.muscle_pass_stiff_scale = ;
% S.subject.tendon_stiff_scale      = ;
S.subject.tendon_stiff_scale = {{'soleus_l','soleus_r','gastroc_r','gastroc_l'},0.7};
% S.subject.mtp_type          = '2022paper';
% S.subject.scale_MT_params         = {{'soleus_l','soleus_r','gastroc_r','gastroc_l'},'FMo',1.2,...
%     {'soleus_l','soleus_r'},'lTs',0.96, {'gastroc_r','gastroc_l'},'lTs',0.985};
% S.subject.spasticity        = ;
% S.subject.muscle_coordination = ;
% S.subject.damping_coefficient_all_dofs = 0; 
S.subject.set_stiffness_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},25};
S.subject.set_damping_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},2};
% S.subject.set_stiffness_coefficient_selected_dofs = {'lumbar_extension',1};
S.subject.set_damping_coefficient_selected_dofs = {{'arm_flex_r','arm_flex_l'},0.2};

S.subject.set_limit_torque_coefficients_selected_dofs = ...
    {{'arm_flex_r','arm_flex_l'},[-10; 22; 10; -22], [-1, 1]};
% S.subject.set_limit_torque_coefficients_selected_dofs = {'all',[0,0,0,0],[0,0]};

% S.subject.base_joints_legs = 'hip';
% S.subject.base_joints_arms = [];

% % S.weights
S.weights.E         = 0.05;
% S.weights.E_exp     = ;
S.weights.q_dotdot  = 1;
S.weights.e_arm     = 10;
S.weights.pass_torq = 0;
S.weights.a         = 1;
S.weights.slack_ctrl = 0.001;
% S.weights.pass_torq_includes_damping = ;

% %S.OpenSimADOptions: required inputs to convert .osim to .dll
% S.OpenSimADOptions.compiler = 'Visual Studio 14 2015 Win64';
S.OpenSimADOptions.verbose_mode = 0; % 0 for no outputs from cmake

        
%% Run predictive simulations

% % warning wrt pelvis heigt for IG
% if S.subject.adapt_IG_pelvis_y == 0 && S.subject.IG_selection ~= "quasi-random"
%     uiwait(msgbox(["Pelvis height of the IG will not be changed.";"Set S.subject.adapt_IG_pelvis_y to 1 if you want to use the model's pelvis height."],"Warning","warn"));
% end

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

