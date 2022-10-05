function [] = PredSim_wrapper_for_app(U,sf)
% inputs
%   * U.ModelName


%% set paths
[pathApp,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathApp);
cd(pathRepo)

%% adapt ModelName
% prevent errors with filenames and paths
ModelName = U.ModelName;
ModelName(strfind(ModelName,' ')) = '_';
ModelName(strfind(ModelName,'/')) = '_';
ModelName(strfind(ModelName,'\')) = '_';

while strcmp(ModelName(1),'_') && length(ModelName)>1
    ModelName = ModelName(2:end);
end

while strcmp(ModelName(end),'_') && length(ModelName)>1
    ModelName = ModelName(1:end-1);
end
U.ModelName = ModelName;

GroupName = U.GroupName;
GroupName(strfind(GroupName,' ')) = '_';
GroupName(strfind(GroupName,'/')) = '_';
GroupName(strfind(GroupName,'\')) = '_';

while strcmp(GroupName(1),'_') && length(GroupName)>1
    GroupName = GroupName(2:end);
end

while strcmp(GroupName(end),'_') && length(GroupName)>1
    GroupName = GroupName(1:end-1);
end
U.GroupName = GroupName;

%% scale model
% Detect height in cm and conver to m
if U.Height > 50
    U.Height = U.Height/100;
end

% Assume constand BMI
U.Mass = 23.1481*U.Height^2;

% Scale
scaleOsim(pathRepo, U, sf)

%% add contact spheres



%% Initialize S
pathDefaultSettings = [pathRepo '\DefaultSettings'];
addpath(pathDefaultSettings)

[S] = initializeSettings();
S.misc.main_path = pathRepo;

addpath([S.misc.main_path '\VariousFunctions'])

%% Required inputs
% name of the subject
S.subject.name = U.ModelName;

% path to folder where you want to store the results of the OCP
S.subject.save_folder  = fullfile('C:\Users\u0150099\OneDrive - KU Leuven\Resultaten_KinderUniversiteit',U.GroupName); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.subject.IG_selection = 'quasi-random';
% S.subject.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Default.mot');
% S.subject.IG_selection_gaitCyclePercent = 50;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = 0;

% casadi folder
S.solver.CasADi_path    = 'C:\GBW_MyPrograms\casadi_3_5_5';


S.subject.v_pelvis_x_trgt   = 1.33;
% S.subject.tendon_stiff      = ;
S.subject.mtp_type          = '2022paper';
S.subject.scale_MT_params         = {{'soleus_l'},'FMo',0.9,};

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
% optional: if you want to install the opensimExe
S.Cpp2Dll.PathCpp2Dll_Exe = InstallOsim2Dll_Exe('C:\GBW_MyPrograms\Osim2Dll_exe'); %(optional: if you want to install the opensimExe)
% S.Cpp2Dll.compiler = 'Visual Studio 15 2017 Win64';
S.Cpp2Dll.export3DSegmentOrigins = [];
S.Cpp2Dll.verbose_mode = 0; % 0 for no outputs from cmake
% S.Cpp2Dll.jointsOrder = ;
% S.Cpp2Dll.coordinatesOrder = ;
        
%% Run predictive simulations
if S.solver.run_as_batch_job
    add_pred_sim_to_batch(S,osim_path)
else
    run_pred_sim(S,osim_path);
end


