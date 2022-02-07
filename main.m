%% Predcitive Simulations of Human Movement

% This script starts the predictive simulation of human movement. The
% required inputs are necessary to start the simulations. Optional inputs,
% if left empty, will be grabed from getDefaultSettings.m.

clear all
close all
clc

%% Initialize S

[S] = initializeSettings();

%% Required inputs

% name of the subject
S.subject.name          = ""; 

% path to folder where you want to store the results of the OCP
S.subject.save_results  = ""; 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.subject.IG_selection  = "";

% give the path to a .mot file on which IG_bounds will be based
S.subject.IG_bounds     = "";

% give the path to the osim model of your subject
osim_model              = "";


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
% S.misc.visualize_IG_bounds = ;
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
% S.solver.max_iter       = ;
% S.solver.parallel_mode  = '';
% S.solver.N_threads      = ;
% S.solver.N_meshes       = ;
% 
% % S.subject
% S.subject.mass              = ;
% s.subject.IG_pelvis_y       = ;
% S.subject.v_pelvis_x_trgt   = ;
% S.subject.muscle_strength   = ;
% S.subject.muscle_stiff      = ;
% S.subject.muscle_sym        = ;
% S.subject.tendon_stiff      = ;
% S.subject.mtp_type          = '';
% S.subject.MT_params         = ;
% S.subject.spasticity        = ;
% S.subject.muscle_coordination = ;
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

