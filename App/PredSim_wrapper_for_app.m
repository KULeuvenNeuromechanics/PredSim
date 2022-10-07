function [] = PredSim_wrapper_for_app(U,sf)
% inputs
%   * U.ModelName


%% set paths
[pathApp,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathApp);
cd(pathRepo);

%% scale model
% Detect height in cm and conver to m
if U.Height > 50
    U.Height = U.Height/100;
end

% Assume constand BMI
% U.Mass = 23.1481*U.Height^2;

% Scale
scaleOsim(pathRepo, U, sf);


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
S.subject.save_folder  = fullfile(U.savefolder,U.GroupName); 

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
% S.subject.IG_selection = 'quasi-random';
S.subject.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Default.mot');
S.subject.IG_selection_gaitCyclePercent = 50;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = 0;

% casadi folder
S.solver.CasADi_path    = 'C:\GBW_MyPrograms\casadi_3_5_5';


S.subject.v_pelvis_x_trgt   = 1.33;
% S.subject.mtp_type          = '2022paper';

S.subject.MT_params  = {{'hamstrings_r','bifemsh_r','glut_max_r','iliopsoas_r',...
    'rect_fem_r','vasti_r','gastroc_r','soleus_r','tib_ant_r','hamstrings_l','bifemsh_l',...
    'glut_max_l','iliopsoas_l','rect_fem_l','vasti_l','gastroc_l','soleus_l','tib_ant_l'},...
    'FMo',U.Force_sf};

S.subject.tendon_stiff_scale      = {{'soleus_l','soleus_r','gastroc_r','gastroc_l'},0.5};

S.subject.set_stiffness_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},25};
S.subject.set_damping_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},2,...
    {'arm_flex_r','arm_flex_l','elbow_flex_r','elbow_flex_l'},0.5};
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
run_pred_sim(S,osim_path);


end