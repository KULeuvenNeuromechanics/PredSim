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
%osimdir = 'C:\OpenSim 4.4\bin\';
%setenv('PATH', [osimdir ';' getenv('PATH')]);


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
S.subject.name = 'Falisse_et_al_2022';

% path to folder where you want to store the results of the OCP
S.subject.save_folder  = fullfile(pathRepoFolder,'PredSimResults',S.subject.name); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
%S.subject.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Full_GC.mot');
S.subject.IG_selection = 'C:\Users\u0130218\Documents\MATLAB\PredSim\Subjects\Falisse_et_al_2022\AverageKineamtics_S204_R.mot';%

%S.subject.IG_selection = "quasi-random";
S.subject.IG_selection_gaitCyclePercent = 100;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% path to folder with program to create dll files from opensim model (will
% be downloaded automatically if it is not there)
%S.Cpp2Dll.PathCpp2Dll_Exe = 'C:\GBW_MyPrograms\Osim2Dll_exe';
S.Cpp2Dll.PathCpp2Dll_Exe = 'C:\Users\u0130218\Documents\MATLAB\osimModelToExtFunc';

% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = 0;

% %%
% % THIS CAN BE USED  IF THE OCP DOESNT OUTPUT THE RESULTS 
% % HAD TO MAKE ADAPTATIONS IN OCP_FORMULATION AS i ADDED 
% % MORE ACTUATORS TO THE KNEE WHICH WERE THOGUHT TO BE OF THE ARM 

% S.post_process.result_filename= 'Falisse_et_al_2022_job1';
% S.post_process.load_prev_opti_vars =1;

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
% S.bounds.coordinates        = {'knee_rot_r',-20,20,'knee_add_r',-20,20, 'knee_rot_l',-20,20,'knee_add_l',-20,20};


% % S.metabolicE - metabolic energy
% S.metabolicE.tanh_b = ;
% S.metabolicE.model  = '';

% % S.misc - miscellanious
% S.misc.v_max_s             = ;
% S.misc.visualize_bounds    = 1;
% S.misc.gaitmotion_type     = '';
% S.misc.msk_geom_eq         = '';
%S.misc.poly_order.lower    = 1;
S.misc.poly_order.upper    = 10 ;
S.misc.msk_geom_bounds      = {{'knee_angle_r','knee_angle_l'},-120,10,{'hip_flexion_r','hip_flexion_l'},-90,90, ...
                                {'ankle_angle_r','ankle_angle_l'},-30,30, {'knee_ty_r','knee_ty_l'},-0.01,0.01 };
% S.misc.gaitmotion_type = ;
% ADED BY BAK 
S.misc.jcf_to_min             = {'knee_ty_r','knee_ty_l'} ;


% % S.post_process
S.post_process.make_plot = 1;
S.post_process.savename  = 'datetime';
% S.post_process.rerun   = 1;
% S.post_process.result_filename = '';

% % S.solver
% S.solver.linear_solver  = '';
% S.solver.tol_ipopt      = ;
% S.solver.max_iter       = 5;
% S.solver.parallel_mode  = '';
% S.solver.N_threads      = 6;
% S.solver.N_meshes       = ;
% S.solver.par_cluster_name = ;
%S.solver.CasADi_path    = 'C:\GBW_MyPrograms\casadi_3_5_5';
S.solver.CasADi_path    = 'C:\Users\u0130218\Documents\MATLAB\Casadi';


% % S.subject
 S.subject.mass              = 75;
% S.subject.IG_pelvis_y       = ;
S.subject.v_pelvis_x_trgt   = 1.33;
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

%S.Cpp2Dll.compiler = 'Visual Studio 17 2022';
%S.Cpp2Dll.compiler = 'Microsoft Visual Studio 12.0';

% S.Cpp2Dll.export3DSegmentOrigins = ;
S.Cpp2Dll.verbose_mode = 0; % 0 for no outputs from cmake
% S.Cpp2Dll.jointsOrder = ;
% S.Cpp2Dll.coordinatesOrder = ;
        
%% Run predictive simulations
% Check for updates in osim2dll
S.Cpp2Dll.PathCpp2Dll_Exe = InstallOsim2Dll_Exe(S.Cpp2Dll.PathCpp2Dll_Exe);


[R] = run_pred_sim_runWIthInitialGuess(S,osim_path);

% 
% % load the reference 
% predSimRef = load('C:\Users\u0130218\Documents\MATLAB\PredSim\Subjects\Falisse_et_al_2022_cleanCopy\InitGuessSim.mat');
% 
% 
% % 
% figure
% plot(R.ground_reaction.GRF_r,'DisplayName','R')
% figure
% plot(predSimRef.R.ground_reaction.GRF_r,'DisplayName','predSimRef')
% 
% % Knee angle R
% figure
% plot(R.kinetics.T_ID(:,10),'r-', 'DisplayName','R')
% hold on
% plot(predSimRef.R.kinetics.T_ID(:,10),'b-',  'DisplayName','predSimRef')
% 
% % Knee angle L
% figure
% plot(R.kinetics.T_ID(:,18),'r-','DisplayName','R')
% hold on
% plot(predSimRef.R.kinetics.T_ID(:,17),'b-','DisplayName','predSimRef')
% 
% % muscel length
% for mi = 14:92
%     h = figure;
%     plot(R.muscles.lMT(:,mi),'r-','DisplayName','R')
%     title(R.model_info.muscle_info.muscle_names(mi))
%     hold on
%     plot(predSimRef.R.muscles.lMT(:,mi),'b-','DisplayName','predSimRef')
%     waitfor(h)
% end
% 
% % ligament lengths
% for li = 1:110
%     h = figure;
%     %plot(R.ligaments.lLig(:,li),'r-','DisplayName','R')
%     plot(R.ligaments.fLig(:,li),'r-','DisplayName','R')
% 
%     title(R.model_info.ligament_info.ligament_names(li))
%  waitfor(h)
% end


