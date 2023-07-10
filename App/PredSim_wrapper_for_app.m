function [varargout] = PredSim_wrapper_for_app(U,sf)
% --------------------------------------------------------------------------
% PredSim_wrapper_for_app
%   This functions replaces main.m when using the Vitruvian_Man apps
% 
% INPUT:
%   - U -
%   * struct with fields
%       U.savefolder    path to high-level folder to save results
%       U.GroupName     name of subfolder to save results
%       U.ModelName     name of the model, 
%       U.Height        heigt of the model (in m)
%       U.Mass          mass of the model (in kg)
%       U.Speed         target speed, or -1 to maximise speed
%
%   - sf -
%   * struct with scale factors
% 
%
% OUTPUT:
%   - (no outputs) -
% 
% Original author: Lars D'Hondt
% Original date: October/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


%% set paths
[pathApp,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathApp);
cd(pathRepo);

%% scale model
% Detect height in cm and convert to m
if U.Height > 50
    U.Height = U.Height/100;
end

% Assume constand BMI
% U.Mass = 23.1481*U.Height^2;

sf_height = U.Height/1.80;
sf_mass = U.Mass/75;

% sf_force = U.Force_sf;
sf_force = sf_mass^(2/3);

% Scale
scaleOsim(pathRepo, U, sf);



%% Inputs
pathDefaultSettings = [pathRepo '\DefaultSettings'];
addpath(pathDefaultSettings)

[S] = initializeSettings();
S.misc.main_path = pathRepo;

addpath([S.misc.main_path '\VariousFunctions'])
addpath([S.misc.main_path '\AdaptOpenSimModel'])

% name of the subject
S.subject.name = U.ModelName;

% path to folder where you want to store the results of the OCP
S.subject.save_folder  = fullfile(U.savefolder,U.GroupName); 

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = 0;

% casadi folder
S.solver.CasADi_path    = U.PathCasadi;
% S.solver.tol_ipopt      = 3;


% scale all with body size (sf_force) and also scale each with segment size^2
S.subject.MT_params  = {
    {'glut_max_r','iliopsoas_r','glut_max_l','iliopsoas_l'},'FMo',sf_force*(sqrt(sf.torso*sf.upp_leg)*sf.upp_leg),...
    {'hamstrings_r','bifemsh_r','rect_fem_r','vasti_r','hamstrings_l','bifemsh_l','rect_fem_l','vasti_l'},'FMo',sf_force*(sf.upp_leg^2),...
    {'gastroc_r','soleus_r','tib_ant_r','gastroc_l','soleus_l','tib_ant_l'},'FMo',sf_force*(sf.low_leg^2)};

% scale all with body size (sf_force) and also scale each with segment size^3 (force * moment arm)
S.subject.scale_actuator_torque = {
    'lumbar_extension',sf_force*(sqrt(sf.torso*sf.shoulder)*sf.torso*sf.upp_leg),...
    {'arm_flex_r','arm_flex_l'},sf_force*(sqrt(sf.torso*sf.shoulder)*sf.shoulder^2),...
    {'elbow_flex_r','elbow_flex_l'},sf_force*sf.upp_arm^3};


S.subject.set_stiffness_coefficient_selected_dofs = {
    {'mtp_angle_l','mtp_angle_r'},25*sf_force*(sf.low_leg^2*sf.foot)};

S.subject.set_damping_coefficient_selected_dofs = {
%     'lumbar_extension',2*sf_force*(sqrt(sf.torso*sf.shoulder)*sf.torso*sf.upp_leg),...
    {'arm_flex_r','arm_flex_l'},0.5*sf_force*sf_force*(sqrt(sf.torso*sf.shoulder)*sf.shoulder^2),...
    {'elbow_flex_r','elbow_flex_l'},0.5*sf_force*sf.upp_arm^3,...
    {'mtp_angle_l','mtp_angle_r'},2*sf_force*(sf.low_leg^2*sf.foot)};


S.subject.tendon_stiff_scale = {
    {'soleus_l','soleus_r','gastroc_r','gastroc_l'},0.7};

S.subject.adapt_IG_pelvis_y = 1;

% % S.weights
% S.weights.E         = 0.05*0;
% % S.weights.E_exp     = ;
% S.weights.q_dotdot  = 1;
% S.weights.e_arm     = 10;
% S.weights.pass_torq = 0;
% S.weights.a         = 1;
% S.weights.slack_ctrl = 0.001;
% S.weights.pass_torq_includes_damping = ;


if U.Speed > 0
    S.subject.v_pelvis_x_trgt = U.Speed;
    S.weights.velocity = 0;
    v_ig = [10:2:20,25:5:45];
    [~,idxv] = min((v_ig/10 - S.subject.v_pelvis_x_trgt).^2);
    name_ig = ['IG_v' num2str(v_ig(idxv)) 'ms.mot'];
    S.subject.IG_selection = fullfile(S.misc.main_path,'Subjects','Vitruvian_Man','IG',name_ig);
    S.subject.IG_selection_gaitCyclePercent = 200;

else
    S.weights.velocity = -5e4;
    S.subject.v_pelvis_x_trgt   = [2,10];

    S.subject.IG_selection = fullfile(S.misc.main_path,'Subjects','Vitruvian_Man','IG','IG_v45ms.mot');
    S.subject.IG_selection_gaitCyclePercent = 200;
end

% S.bounds.coordinates = [];

S.OpenSimADOptions.verbose_mode = false;
S.OpenSimADOptions.compiler = 'Visual Studio 14 2015 Win64';

%% Run predictive simulations
savename = run_pred_sim(S,osim_path);

savepath = fullfile(S.subject.save_folder,savename);

if nargout>=1
    varargout{1} = savepath;
end

end