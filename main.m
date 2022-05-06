%% Predcitive Simulations of Human Movement

% This script starts the predictive simulation of human movement. The
% required inputs are necessary to start the simulations. Optional inputs,
% if left empty, will be taken from getDefaultSettings.m.

clear all
close all
clc
[pathRepo,~,~] = fileparts(mfilename('fullpath'));

%% Initialize S
pathDefaultSettings = [pathRepo '\DefaultSettings'];
addpath(pathDefaultSettings)

[S] = initializeSettings();
S.misc.main_path = pathRepo;

%% Required inputs
% name of the subject
S.subject.name          = 'test_1'; 

% path to folder where you want to store the results of the OCP
S.subject.save_folder  = fullfile(pathRepo,'test_1_results'); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
% S.subject.IG_selection  = 'C:\Users\u0150099\Documents\master_thesis\ReferenceData\ModelScaling\reference_data\IK_subject1_withMTJ_locked_scaled_default_gait_14.mot';
% S.subject.IG_selection_gaitCyclePercent = 100;
S.subject.IG_selection = 'quasi-random';

% give the path to the osim model of your subject
osim_path              = fullfile(pathRepo,'Subjects','test_1','test_1.osim');


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
% 
% % S.metabolicE - metabolic energy
% S.metabolicE.tanh_b = ;
% S.metabolicE.model  = '';
% 
% % S.misc - miscellanious
% S.misc.v_max_s             = ;
% S.misc.visualize_IG_bounds = 1;
% S.misc.gaitmotion_type     = '';
% S.misc.msk_geom_eq         = '';
% S.misc.poly_order.lower    = ;
% S.misc.poly_order.upper    = ;
% 
% % S.post_process
% S.post_process.make_plot = '';
% S.post_process.savename  = '';
% 
% % S.solver
% S.solver.linear_solver  = '';
% S.solver.tol_ipopt      = ;
S.solver.max_iter       = 5;
% S.solver.parallel_mode  = '';
% S.solver.N_threads      = ;
% S.solver.N_meshes       = ;
% 
% % S.subject
S.subject.mass              = 75;
S.subject.IG_pelvis_y       = 0.89;
S.subject.v_pelvis_x_trgt   = 1.33;
% S.subject.IG_bounds = ;
% S.subject.muscle_strength   = ;
% S.subject.muscle_pass_stiff_shift = ;
% S.subject.muscle_pass_stiff_scale = ;
% S.subject.muscle_sym        = 1;
% S.subject.tendon_stiff      = ;
% S.subject.mtp_type          = 'passive';
% S.subject.MT_params         = ;
% S.subject.spasticity        = ;
% S.subject.muscle_coordination = ;
S.subject.set_stiffness_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},25};
S.subject.set_damping_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},2};
% 
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

run_pred_sim(S,osim_path);

%% Plot figures

plot_figures

