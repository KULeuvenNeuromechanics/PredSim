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

%% Required inputs
% name of the subject
S.subject.name = 'CP3_T0_scaled_MRI_v7_scaledMT_test'; 

% path to folder where you want to store the results of the OCP
S.subject.save_folder  = fullfile(pathRepoFolder,'PredSimResults',S.subject.name); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
% S.subject.IG_selection = 'quasi-random';
S.subject.IG_selection = fullfile(S.misc.main_path,'OCP','Fal_s1_mtp_v2.mot');
S.subject.IG_selection_gaitCyclePercent = 200;

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
% S.bounds.coordinates        = {'pelvis_tilt',-30*pi/180,30*pi/180};

% % S.metabolicE - metabolic energy
% S.metabolicE.tanh_b = ;
% S.metabolicE.model  = '';

% % S.misc - miscellanious
% S.misc.v_max_s             = ;
S.misc.visualize_bounds = 1;
% S.misc.gaitmotion_type     = '';
% S.misc.msk_geom_eq         = '';
% S.misc.poly_order.lower    = ;
% S.misc.poly_order.upper    = ;

% % S.post_process
% S.post_process.make_plot = '';
% S.post_process.savename  = '';
% S.post_process.rerun   = 1;
% S.post_process.result_filename = 'Fal_s1_mtp_v2';

% % S.solver
% S.solver.linear_solver  = '';
% S.solver.tol_ipopt      = ;
% S.solver.max_iter       = 5;
% S.solver.parallel_mode  = '';
S.solver.N_threads      = 8;
% S.solver.N_meshes       = ;
% S.solver.par_cluster_name = ;
S.solver.CasADi_path    = 'C:\GBW_MyDownloads\CasADi v3_5_2';


% % S.subject
S.subject.mass              = 30.15;
S.subject.IG_pelvis_y       = 0.9385;
S.subject.v_pelvis_x_trgt   = 1.33;
% S.subject.IG_bounds = ;
% S.subject.muscle_strength   = ;
% S.subject.muscle_pass_stiff_shift = ;
% S.subject.muscle_pass_stiff_scale = ;
% S.subject.muscle_sym        = 1;
% S.subject.tendon_stiff      = ;
S.subject.mtp_type          = '2022paper';
% S.subject.MT_params         = ;
% S.subject.spasticity        = ;
% S.subject.muscle_coordination = ;
S.subject.set_stiffness_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},17};
S.subject.set_damping_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},0.5};
S.subject.set_limit_torque_coefficients_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},[0,0,0,0],[0,0]};

% % S.weights
% S.weights.E         = ;
% S.weights.E_exp     = ;
% S.weights.q_dotdot  = ;
% S.weights.e_arm     = ;
% S.weights.pass_torq = ;
% S.weights.a         = ;
% S.weights.e_mtp     = ;
% S.weights.slack_ctrl = ;


%% Run predictive simulations
if S.solver.run_as_batch_job
    add_pred_sim_to_batch(S,osim_path)
else
    run_pred_sim(S,osim_path);
end


%% Plot figures

% plot_figures

