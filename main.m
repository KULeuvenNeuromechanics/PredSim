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

[S] = initializeSettings('Falisse_et_al_2022');
S.misc.main_path = pathRepo;

addpath(fullfile(S.misc.main_path,'VariousFunctions'))

%% Required inputs
% name of the subject
S.subject.name = 'Falisse_et_al_2022';

% path to folder where you want to store the results of the OCP
S.subject.save_folder  = fullfile(pathRepoFolder,'PredSimResults',S.subject.name); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
% S.subject.IG_selection = 'quasi-random';
S.subject.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
S.subject.IG_selection_gaitCyclePercent = 100;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = 0;

%% Optional inputs
% see README.md in the main folder for information about these optional
% inputs.

% % S.bounds
% S.bounds.a.lower            = ;
% S.bounds.SLL.upper          = ;
% S.bounds.SLR.upper          = ;
% S.bounds.dist_trav.lower    = ;
% S.bounds.t_final.upper      = ;
% S.bounds.t_final.lower      = ;
% S.bounds.Qs                 = {'pelvis_tilt',-30,30,'pelvis_list',-30,30};



% % S.metabolicE - metabolic energy
% S.metabolicE.tanh_b = 100;
% S.metabolicE.model  = '';

% % S.misc - miscellanious
% S.misc.v_max_s             = ;
% S.misc.visualize_bounds    = 1;
% S.misc.gaitmotion_type     = '';
% S.misc.msk_geom_eq         = '';
% S.misc.poly_order.lower    = ;
% S.misc.poly_order.upper    = ;
% S.misc.default_msk_geom_bound = ;
% S.misc.msk_geom_bounds      = {{'knee_angle_r','knee_angle_l'},-120,10,'lumbar_extension',nan,30};
% S.misc.gaitmotion_type = 'FullGaitCycle';

% % S.post_process
S.post_process.make_plot = 1;
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
S.solver.N_meshes       = 50;
% S.solver.par_cluster_name = ;
S.solver.CasADi_path = casadi.GlobalOptions.getCasadiPath(); % ask casadi



% % S.subject
% S.subject.mass              = ;
% S.subject.IG_pelvis_y       = 1;
% S.subject.adapt_IG_pelvis_y = 1;
% S.subject.v_pelvis_x_trgt   = 1.33;
% S.subject.muscle_strength   = ;
% S.subject.muscle_pass_stiff_shift = {{'soleus','_gas','per_','tib_','_dig_','_hal_','FDB'},0.9}; %,'FDB'
% S.subject.muscle_pass_stiff_scale = ;
% S.subject.tendon_stiff_scale      = {{'soleus','gastroc'},0.5};
% S.subject.scale_MT_params = {{'soleus_l'},'FMo',0.9,{'soleus_l'},'alphao',1.1};
% S.subject.spasticity        = ;
% S.subject.muscle_coordination = ;
% S.subject.mtp_type          = '2022paper';
% S.subject.set_stiffness_coefficient_selected_dofs = {'mtp_angle',1};
% S.subject.set_damping_coefficient_selected_dofs = {'mtp_angle',2};
% S.subject.set_limit_torque_coefficients_selected_dofs = {{'mtj_angle_l','mtj_angle_r'},[0,0,0,0],[0,0]};
% S.subject.base_joints_arms = [];


% % S.weights
% S.weights.E         = 0;
% % S.weights.E_exp     = ;
% S.weights.q_dotdot  = 1;
% S.weights.e_arm     = 500;
% S.weights.pass_torq = 1;
% S.weights.a         = 10*18;
% S.weights.slack_ctrl = ;


% % S.orthosis
% add orthosis on right side
S.orthosis.settings{1}.function_name = 'parametricAFO';
S.orthosis.settings{1}.ankle_stiffness = 100; % Nm/rad
S.orthosis.settings{1}.mtp_stiffness = 10; % Nm/rad
S.orthosis.settings{1}.left_right = 'r';

% add the same orthosis on left side
S.orthosis.settings{2}.function_name = 'parametricAFO';
S.orthosis.settings{2}.ankle_stiffness = 100; % Nm/rad
S.orthosis.settings{2}.mtp_stiffness = 10; % Nm/rad
S.orthosis.settings{2}.left_right = 'l';


% %S.OpenSimADOptions: required inputs to convert .osim to .dll
% S.OpenSimADOptions.compiler = 'Visual Studio 17 2022';
S.OpenSimADOptions.verbose_mode = 0; % 0 for no outputs from cmake


%% Run predictive simulations

% warning wrt pelvis heigt for IG
if S.subject.adapt_IG_pelvis_y == 0 && S.subject.IG_selection ~= "quasi-random"
%     uiwait(msgbox(["Pelvis height of the IG will not be changed.";"Set S.subject.adapt_IG_pelvis_y to 1 if you want to use the model's pelvis height."],"Warning","warn"));
end

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

