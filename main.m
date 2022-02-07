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

Main % user specified settings and osim path
    []=run_pred_sim(S,osim_path)
        S=getDefaultSettings(S) % adding settings that are not specified
        %%
        [S,model_info]=preprocessing(S,osim_path)
            []=osim2dll(S,osim_path) % save ['F_' osim_file_name '.dll'] and ['F_' osim_file_name '_IO.mat'] in subject-folder (.mat in 32bit int)
            model_info = get_model_info(S,osim_path) % load IO and create struct model_info (convert 32bit to double)
            [model_info] = read_and_scale_MTparameters(S,osim_path,model_info) % save MTparameters in subject-folder
                MTparameters = getMTparameters(osim_path,muscle_names) % grab muscle_names from model_info
                *MTparameters = scale_MTparameters(S,MTparameters) % scale lM0,vM0,lTs,FMopt,alpha0 (as specified in S)
            % only run if the geometry for the given model does not exist yet
            get_musculoskeletal_geometry_approximation(S,osim_path,model_info)
                muscle_data = muscle_analysis(S,osim_path,model_info) 
                [muscle_spanning_joint_info, muscle_info] = PolynomialFit(muscle_data) 
                % save muscle_spanning_joint_INFO, muscle_info, muscle_data in subject-folder
            model_info = update_model_info(model_info) % https://github.com/KULeuvenNeuromechanics/PredSim/blob/workProgrammingRetreat/OCP/Copy_of_f_PredSim_Gait92.m#L34
        %%
        [f_casadi] = createCasadiFunctions(main_path,model_info) % number of inputs will changed based on pre-processing
        % https://github.com/KULeuvenNeuromechanics/PredSim/blob/dev-buurke/CasadiFunctions/createCasadiFunctions.m#L1
            % need to check the functions that are already in this folder
        %%
        []=OCP_formulation(S,model_info,f_casadi) 
        % https://github.com/KULeuvenNeuromechanics/PredSim/blob/workProgrammingRetreat/OCP/Copy_of_f_PredSim_Gait92.m
        % need to finish based on pre-processing and casadifunctions
            % split into 2 parts:
                % 1) common for formulation and post-processing -> run every time
                % 2) make collocation casadi formulation -> run only if solving
            % extract results, assert cost, reconstruct gait cycle
            % save results in specified folder
        %%
        []=post_processing(S,model_info,f_casadi)
            % load results from specified folder
            % specific things
                % e.g.: [R] = getGRFS(R,S,modelinfo) to extract the GRFs from the external function
                % save R
        %%
plot_figures % plot info for a list of R

% next meeting: validation steps (early january)

* not crucial to run


%% notes

- validation flag for each function: is it validated or not















