function [] = f_PredSim_PoggenSee2020_DefaultS(S)


%% Notes
rmpath(genpath('G:\PhD\Matlabcodes\Casadi\casadi'));
addpath(genpath('D:\MaartenAfschrift\softInstall\casadi'));
addpath(genpath('D:\MaartenAfschrift\GitProjects\3dpredictsim'));

% to simplify batch processing, the casadi functions were already created
% using the script CasadiFunctions_all_mtp_createDefault.m
% This assumes invariant:
%   Muscle and tendon properties
%   Polynomials to compute moment arms
%   Functions to compute passive stiffness

% We can still vary:
% 1) the collocation scheme
% 2) the weights in the objective function
% 3) the exoskeleton assistance
% 4) the external function

%% Default settings

if ~isfield(S,'savename_ig')
    S.savename_ig= [];
end

if ~isfield(S,'ResultsF_ig')
    S.ResultsF_ig = [];
end

if  ~isfield(S,'NThreads') || isempty(S.NThreads)
    S.NThreads = 4;
end

if ~isfield(S,'mass') || isempty(S.mass)
    S.mass = 62;
end

if ~isfield(S,'subject') || isempty(S.subject)
    S.subject = 'subject1';
end

% quasi random initial guess
if ~isfield(S,'IG_PelvisY') || isempty(S.IG_PelvisY)
    S.IG_PelvisY = 0.9385;
end

%% User inputs (typical settings structure)
% load default CasadiFunctions

% flow control
solveProblem    = S.Flow.solveProblem; % set to 1 to solve problem

% settings for optimization
v_tgt       = S.v_tgt;      % average speed
N           = S.N;          % number of mesh intervals
W.E         = S.W.E;        % weight metabolic energy rate
W.Ak        = S.W.Ak;       % weight joint accelerations
W.ArmE      = S.W.ArmE;     % weight arm excitations
W.passMom   = S.W.passMom;  % weight passive torques
W.A         = S.W.A;        % weight muscle activations
exp_E       = S.W.exp_E;    % power metabolic energy
W.Mtp       = S.W.Mtp;      % weight mtp excitations
W.u         = S.W.u;        % weight on exctiations arm actuators
IGsel       = S.IGsel;      % initial guess identifier
IGm         = S.IGmodeID;   % initial guess mode identifier
coCont      = S.coCont;     % co-contraction identifier

% identifier for EMG load
savename_ig = S.savename_ig;

% ipopt options
tol_ipopt       = S.tol_ipopt;
linear_solver   = S.linear_solver;

%% Settings

import casadi.*
parallelMode = 'thread';
NThreads = S.NThreads;
subject = S.subject;


%% Load external functions
% The external function performs inverse dynamics through the
% OpenSim/Simbody C++ API. This external function is compiled as a dll from
% which we create a Function instance using CasADi in MATLAB. More details
% about the external function can be found in the documentation.
pathmain = pwd;
% We use different external functions, since we also want to access some
% parameters of the model in a post-processing phase.
[pathRepo,~,~] = fileparts(pathmain);
pathExternalFunctions = [pathRepo,'/ExternalFunctions'];
% Loading external functions.
cd(pathExternalFunctions);
setup.derivatives =  'AD'; % Algorithmic differentiation
F  = external('F',S.ExternalFunc);
cd(pathmain);

%% Output folder
% save the results in the right folder:
OutFolder = fullfile(pathRepo,'Results',S.ResultsFolder);
if ~isfolder(OutFolder)
    mkdir(OutFolder);
end

%% Indices external function
% Indices of the elements in the external functions
% External function: F
% First, joint torques.
jointi.pelvis.tilt  = 1;
jointi.pelvis.list  = 2;
jointi.pelvis.rot   = 3;
jointi.pelvis.tx    = 4;
jointi.pelvis.ty    = 5;
jointi.pelvis.tz    = 6;
jointi.hip_flex.l   = 7;
jointi.hip_add.l    = 8;
jointi.hip_rot.l    = 9;
jointi.hip_flex.r   = 10;
jointi.hip_add.r    = 11;
jointi.hip_rot.r    = 12;
jointi.knee.l       = 13;
jointi.knee.r       = 14;
jointi.ankle.l      = 15;
jointi.ankle.r      = 16;
jointi.subt.l       = 17;
jointi.subt.r       = 18;
jointi.mtp.l        = 19;
jointi.mtp.r        = 20;
jointi.trunk.ext    = 21;
jointi.trunk.ben    = 22;
jointi.trunk.rot    = 23;
jointi.sh_flex.l    = 24;
jointi.sh_add.l     = 25;
jointi.sh_rot.l     = 26;
jointi.sh_flex.r    = 27;
jointi.sh_add.r     = 28;
jointi.sh_rot.r     = 29;
jointi.elb.l        = 30;
jointi.elb.r        = 31;
% Vectors of indices for later use
residualsi          = jointi.pelvis.tilt:jointi.elb.r; % all
ground_pelvisi      = jointi.pelvis.tilt:jointi.pelvis.tz; % ground-pelvis
trunki              = jointi.trunk.ext:jointi.trunk.rot; % trunk
armsi               = jointi.sh_flex.l:jointi.elb.r; % arms
mtpi                = jointi.mtp.l:jointi.mtp.r; % mtps
residuals_noarmsi   = jointi.pelvis.tilt:jointi.trunk.rot; % all but arms
roti                = [jointi.pelvis.tilt:jointi.pelvis.rot,...
    jointi.hip_flex.l:jointi.elb.r];
% Number of degrees of freedom for later use
nq.all      = length(residualsi); % all
nq.abs      = length(ground_pelvisi); % ground-pelvis
nq.trunk    = length(trunki); % trunk
nq.arms     = length(armsi); % arms
nq.mtp     = length(mtpi); % arms
nq.leg      = 10; % #joints needed for polynomials
% Second, origins bodies.
% Calcaneus
calcOr.r    = 32:33;
calcOr.l    = 34:35;
calcOr.all  = [calcOr.r,calcOr.l];
NcalcOr     = length(calcOr.all);
% Femurs
femurOr.r   = 36:37;
femurOr.l   = 38:39;
femurOr.all = [femurOr.r,femurOr.l];
NfemurOr    = length(femurOr.all);
% Hands
handOr.r    = 40:41;
handOr.l    = 42:43;
handOr.all  = [handOr.r,handOr.l];
NhandOr     = length(handOr.all);
% Tibias
tibiaOr.r   = 44:45;
tibiaOr.l   = 46:47;
tibiaOr.all = [tibiaOr.r,tibiaOr.l];
NtibiaOr    = length(tibiaOr.all);
% External function: F1 (post-processing purpose only)
% Ground reaction forces (GRFs)
GRFi.r      = 32:34;
GRFi.l      = 35:37;
GRFi.all    = [GRFi.r,GRFi.l];
NGRF        = length(GRFi.all);
% Origins calcaneus (3D)
calcOrall.r     = 38:40;
calcOrall.l     = 41:43;
calcOrall.all   = [calcOrall.r,calcOrall.l];
NcalcOrall      = length(calcOrall.all);

%% Model info
body_weight = S.mass*9.81;

%% Collocation scheme
% We use a pseudospectral direct collocation method, i.e. we use Lagrange
% polynomials to approximate the state derivatives at the collocation
% points in each mesh interval. We use d=3 collocation points per mesh
% interval and Radau collocation points.
pathCollocationScheme = [pathRepo,'/CollocationScheme'];
addpath(genpath(pathCollocationScheme));
d = 3; % degree of interpolating polynomial
method = 'radau'; % collocation method
[tau_root,C,D,B] = CollocationScheme(d,method);

%% Muscle-tendon parameters
% Muscles from one leg and from the back
muscleNames = {'glut_med1_r','glut_med2_r','glut_med3_r',...
    'glut_min1_r','glut_min2_r','glut_min3_r','semimem_r',...
    'semiten_r','bifemlh_r','bifemsh_r','sar_r','add_long_r',...
    'add_brev_r','add_mag1_r','add_mag2_r','add_mag3_r','tfl_r',...
    'pect_r','grac_r','glut_max1_r','glut_max2_r','glut_max3_r',......
    'iliacus_r','psoas_r','quad_fem_r','gem_r','peri_r',...
    'rect_fem_r','vas_med_r','vas_int_r','vas_lat_r','med_gas_r',...
    'lat_gas_r','soleus_r','tib_post_r','flex_dig_r','flex_hal_r',...
    'tib_ant_r','per_brev_r','per_long_r','per_tert_r','ext_dig_r',...
    'ext_hal_r','ercspn_r','intobl_r','extobl_r','ercspn_l',...
    'intobl_l','extobl_l'};
% Muscle indices for later use
pathmusclemodel = [pathRepo,'/MuscleModel'];
addpath(genpath(pathmusclemodel));
% (1:end-3), since we do not want to count twice the back muscles
musi = MuscleIndices(muscleNames(1:end-3));
% Total number of muscles
NMuscle = length(muscleNames(1:end-3))*2;
% Muscle-tendon parameters. Row 1: maximal isometric forces; Row 2: optimal
% fiber lengths; Row 3: tendon slack lengths; Row 4: optimal pennation
% angles; Row 5: maximal contraction velocities
pathpolynomial = fullfile(pathRepo,'Polynomials',subject);
addpath(genpath(pathpolynomial));
tl = load([pathpolynomial,'/muscle_spanning_joint_INFO_',subject,'_mtp.mat']);
[~,mai] = MomentArmIndices(muscleNames(1:end-3),...
    tl.muscle_spanning_joint_INFO(1:end-3,:));

% Parameters for activation dynamics
tact = 0.015; % Activation time constant
tdeact = 0.06; % Deactivation time constant

%% Metabolic energy model parameters
% We extract the specific tensions and slow twitch rations.
pathMetabolicEnergy = [pathRepo,'/MetabolicEnergy'];
addpath(genpath(pathMetabolicEnergy));
% (1:end-3), since we do not want to count twice the back muscles
tension = getSpecificTensions(muscleNames(1:end-3));
tensions = [tension;tension];
% (1:end-3), since we do not want to count twice the back muscles
pctst = getSlowTwitchRatios(muscleNames(1:end-3));
pctsts = [pctst;pctst];

%% CasADi functions
% We create several CasADi functions for later use
pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
PathDefaultFunc = fullfile(pathCasADiFunctions,S.CasadiFunc_Folders);
cd(PathDefaultFunc);
% f_coll = Function.load('f_coll');
f_ArmActivationDynamics = Function.load('f_ArmActivationDynamics');
f_FiberLength_TendonForce_tendon = Function.load('f_FiberLength_TendonForce_tendon');
f_FiberVelocity_TendonForce_tendon = Function.load('f_FiberVelocity_TendonForce_tendon');
f_forceEquilibrium_FtildeState_all_tendon = Function.load('f_forceEquilibrium_FtildeState_all_tendon');
f_J2    = Function.load('f_J2');
f_J23   = Function.load('f_J23');
f_J25   = Function.load('f_J25');
f_J8    = Function.load('f_J8');
f_J92   = Function.load('f_J92');
f_J92exp = Function.load('f_J92exp');
f_Jnn2  = Function.load('f_Jnn2');
f_Jnn3  = Function.load('f_Jnn3');
f_lMT_vMT_dM = Function.load('f_lMT_vMT_dM');
f_MtpActivationDynamics = Function.load('f_MtpActivationDynamics');
% f_PassiveMoments = Function.load('f_PassiveMoments');
% f_passiveTATorque = Function.load('f_passiveTATorque');
f_T12 = Function.load('f_T12');
f_T13 = Function.load('f_T13');
f_T27 = Function.load('f_T27');
f_T6 = Function.load('f_T6');
f_AllPassiveTorques = Function.load('f_AllPassiveTorques');
fgetMetabolicEnergySmooth2004all = Function.load('fgetMetabolicEnergySmooth2004all');
cd(pathmain);


%% Experimental data
% We extract experimental data to set bounds and initial guesses if needed
pathData = [pathRepo,'/OpenSimModel/',subject];
joints = {'pelvis_tilt','pelvis_list','pelvis_rotation','pelvis_tx',...
    'pelvis_ty','pelvis_tz','hip_flexion_l','hip_adduction_l',...
    'hip_rotation_l','hip_flexion_r','hip_adduction_r','hip_rotation_r',...
    'knee_angle_l','knee_angle_r','ankle_angle_l','ankle_angle_r',...kbk
    'subtalar_angle_l','subtalar_angle_r','mtp_angle_l','mtp_angle_r',...
    'lumbar_extension','lumbar_bending','lumbar_rotation','arm_flex_l',...
    'arm_add_l','arm_rot_l','arm_flex_r','arm_add_r','arm_rot_r',...
    'elbow_flex_l','elbow_flex_r'};
pathVariousFunctions = [pathRepo,'/VariousFunctions'];
addpath(genpath(pathVariousFunctions));
% Extract joint positions from average walking motion
motion_walk         = 'walking';
nametrial_walk.id   = ['average_',motion_walk,'_HGC_mtp'];
nametrial_walk.IK   = ['IK_',nametrial_walk.id];
pathIK_walk         = [pathData,'/IK/',nametrial_walk.IK,'.mat'];
Qs_walk             = getIK(pathIK_walk,joints);
% Depending on the initial guess mode, we extract further experimental data
if IGm == 2
    % Extract joint positions from average running motion
    motion_run          = 'running';
    nametrial_run.id    = ['average_',motion_run,'_HGC'];
    nametrial_run.IK    = ['IK_',nametrial_run.id];
    pathIK_run          = [pathData,'/IK/',nametrial_run.IK,'.mat'];
    Qs_run              = getIK(pathIK_run,joints);
elseif IGm == 3
    % Extract joint positions from existing motion (previous results)
    GuessFolder = fullfile(pathRepo,'Results',S.ResultsF_ig);
    pathIK      = fullfile(GuessFolder,[savename_ig '.mot']);
    Qs_ig       = getIK(pathIK,joints);
    % When saving the results, we save a full gait cycle (2*N) so here we
    % only select 1:N to have half a gait cycle
    Qs_ig_sel.allfilt   = Qs_ig.allfilt(1:N,:);
    Qs_ig_sel.time      = Qs_ig.time(1:N,:);
    Qs_ig_sel.colheaders = Qs_ig.colheaders;
end

%% Bounds
pathBounds = [pathRepo,'/Bounds'];
addpath(genpath(pathBounds));
[bounds,scaling] = getBounds_all_mtp(Qs_walk,NMuscle,nq,jointi,v_tgt);
% Simulate co-contraction by increasing the lower bound on muscle activations
if coCont == 1
    bounds.a.lower = 0.1*ones(1,NMuscle);
elseif coCont == 2
    bounds.a.lower = 0.15*ones(1,NMuscle);
elseif coCont == 3
    bounds.a.lower = 0.2*ones(1,NMuscle);
end

%% Initial guess
% The initial guess depends on the settings
pathIG = [pathRepo,'/IG'];
addpath(genpath(pathIG));
if IGsel == 1 % Quasi-random initial guess
    guess = getGuess_QR_opti_int(N,nq,NMuscle,scaling,v_tgt,jointi,d,S.IG_PelvisY);
elseif IGsel == 2 % Data-informed initial guess
    if IGm == 1 % Data from average walking motion
        time_IC = [Qs_walk.time(1),Qs_walk.time(end)];
        guess = getGuess_DI_opti_int_mtp(Qs_walk,nq,N,time_IC,NMuscle,jointi,...
            scaling,v_tgt,d);
    elseif IGm == 2 % Data from average runing motion
        time_IC = [Qs_run.time(1),Qs_run.time(end)];
        guess = getGuess_DI_opti_int(Qs_run,nq,N,time_IC,NMuscle,jointi,...
            scaling,v_tgt,d);
    elseif IGm == 3 % Data from selected motion
        time_IC = [Qs_ig_sel.time(1),Qs_ig_sel.time(end)];
        guess = getGuess_DI_opti_int_mtp(Qs_ig_sel,nq,N,time_IC,NMuscle,jointi,scaling,v_tgt,d);
    end
end
% If co-contraction, the initial guess of muscles activations is increased
if coCont == 1
    guess.a = 0.15*ones(N,NMuscle);
elseif coCont == 2
    guess.a = 0.20*ones(N,NMuscle);
elseif coCont == 3
    guess.a = 0.25*ones(N,NMuscle);
end

%% exoskeleton torques
ExoControl = [];
body_mass = S.mass;
if S.ExoBool
    if strcmp(S.DataSet,'Zhang2017')
        % load the data from Zhang 2017
        [DataPath,~,~] = fileparts(pathRepo);
        Zhang = load([DataPath,'\Data\Zhang_2017\opt_tau.mat']);
        Tankle = nanmean(Zhang.opt_tau)*-1.*body_mass; % -1 because plantarflexion is negative in opensim model
        ExoSpline.Tankle = spline(linspace(0,2,length(Tankle)*2),[Tankle Tankle]);
    elseif strcmp(S.DataSet,'PoggenSee2020_AFO')
        [DataPath,~,~] = fileparts(pathRepo);
        Poggensee = load([DataPath,'\Data\Poggensee_2020\torque_profile.mat']);
        Tankle = Poggensee.torque*-1*body_mass; % -1 because plantarflexion is negative in opensim model
        ExoSpline.Tankle = spline(linspace(0,2,length(Tankle)*2),[Tankle' Tankle']);
    else
        error(['Could not find the dataset ' S.DataSet ' to prescribe the exoskeleton torques']);
    end
    
    ExoControl.Tankle_r = ppval(ExoSpline.Tankle,linspace(0,0.5,N));
    ExoControl.Tankle_l = ppval(ExoSpline.Tankle,linspace(0.5,1,N));
    if isfield(S,'ExoScale')
        ExoControl.Tankle_r = ExoControl.Tankle_r*S.ExoScale;
        ExoControl.Tankle_l = ExoControl.Tankle_l*S.ExoScale;
    end
    ExoVect = [ExoControl.Tankle_l; ExoControl.Tankle_r];
else
    ExoVect = zeros(2,N);
end

%% Index helpers

% indexes to select kinematics left and right leg
IndexLeft = [jointi.hip_flex.l jointi.hip_add.l jointi.hip_rot.l, ...
    jointi.knee.l jointi.ankle.l jointi.subt.l jointi.mtp.l,...
    jointi.trunk.ext, jointi.trunk.ben, jointi.trunk.rot];
IndexRight = [jointi.hip_flex.r jointi.hip_add.r jointi.hip_rot.r, ...
    jointi.knee.r jointi.ankle.r jointi.subt.r jointi.mtp.r,...
    jointi.trunk.ext, jointi.trunk.ben, jointi.trunk.rot];

% indexes for symmetry steps
QsInvA = [jointi.pelvis.tilt,...
    jointi.pelvis.ty,...
    jointi.hip_flex.l:jointi.trunk.ext,...
    jointi.sh_flex.l:jointi.elb.r]';
QsInvB = [jointi.pelvis.tilt,...
    jointi.pelvis.ty,...
    jointi.hip_flex.r:jointi.hip_rot.r,...
    jointi.hip_flex.l:jointi.hip_rot.l,...
    jointi.knee.r,...
    jointi.knee.l,...
    jointi.ankle.r,...
    jointi.ankle.l,...
    jointi.subt.r,...
    jointi.subt.l,...
    jointi.mtp.r,...
    jointi.mtp.l,...
    jointi.trunk.ext,...
    jointi.sh_flex.r:jointi.sh_rot.r,...
    jointi.sh_flex.l:jointi.sh_rot.l,...
    jointi.elb.r,...
    jointi.elb.l]';

QdotsInvA = [jointi.pelvis.tilt,...
    jointi.pelvis.tx,jointi.pelvis.ty,...
    jointi.hip_flex.l:jointi.trunk.ext,...
    jointi.sh_flex.l:jointi.elb.r]';
QdotsInvB = [jointi.pelvis.tilt,...
    jointi.pelvis.tx,jointi.pelvis.ty,...
    jointi.hip_flex.r:jointi.hip_rot.r,...
    jointi.hip_flex.l:jointi.hip_rot.l,...
    jointi.knee.r,...
    jointi.knee.l,...
    jointi.ankle.r,...
    jointi.ankle.l,...
    jointi.subt.r,...
    jointi.subt.l,...
    jointi.mtp.r,...
    jointi.mtp.l,...
    jointi.trunk.ext,...
    jointi.sh_flex.r:jointi.sh_rot.r,...
    jointi.sh_flex.l:jointi.sh_rot.l,...
    jointi.elb.r,...
    jointi.elb.l]';

orderQsOpp = [jointi.pelvis.list:jointi.pelvis.list,...
    jointi.pelvis.rot:jointi.pelvis.rot,...
    jointi.pelvis.tz:jointi.pelvis.tz,...
    jointi.trunk.ben:jointi.trunk.ben,...
    jointi.trunk.rot:jointi.trunk.rot];

% Helper variables to reconstruct full gait cycle assuming symmetry
QsSymA = [jointi.pelvis.tilt,jointi.pelvis.ty,...
    jointi.hip_flex.l:jointi.trunk.ext,...
    jointi.sh_flex.l:jointi.elb.r];
QsSymB = [jointi.pelvis.tilt,jointi.pelvis.ty,...
    jointi.hip_flex.r:jointi.hip_rot.r,...
    jointi.hip_flex.l:jointi.hip_rot.l,...
    jointi.knee.r,jointi.knee.l,...
    jointi.ankle.r,jointi.ankle.l,...
    jointi.subt.r,jointi.subt.l,...
    jointi.mtp.r,jointi.mtp.l,...
    jointi.trunk.ext,...
    jointi.sh_flex.r:jointi.sh_rot.r,...
    jointi.sh_flex.l:jointi.sh_rot.l,...
    jointi.elb.r,jointi.elb.l];
QsOpp = [jointi.pelvis.list:jointi.pelvis.rot,jointi.pelvis.tz,...
    jointi.trunk.ben:jointi.trunk.rot];
QsSymA_ptx = [jointi.pelvis.tilt,jointi.pelvis.tx,...
    jointi.pelvis.ty,...
    jointi.hip_flex.l:jointi.trunk.ext,...
    jointi.sh_flex.l:jointi.elb.r];
QsSymB_ptx = [jointi.pelvis.tilt,jointi.pelvis.tx,...
    jointi.pelvis.ty,...
    jointi.hip_flex.r:jointi.hip_rot.r,...
    jointi.hip_flex.l:jointi.hip_rot.l,...
    jointi.knee.r,jointi.knee.l,...
    jointi.ankle.r,jointi.ankle.l,...
    jointi.subt.r,jointi.subt.l,...
    jointi.mtp.r,jointi.mtp.l,...
    jointi.trunk.ext,...
    jointi.sh_flex.r:jointi.sh_rot.r,...
    jointi.sh_flex.l:jointi.sh_rot.l,...
    jointi.elb.r,jointi.elb.l];



if solveProblem
    
    %% OCP create variables and bounds
    
    opti = casadi.Opti();
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define static parameters
    % Final time
    tf = opti.variable();
    opti.subject_to(bounds.tf.lower < tf < bounds.tf.upper);
    opti.set_initial(tf, guess.tf);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define states
    % Muscle activations at mesh points
    a = opti.variable(NMuscle,N+1);
    opti.subject_to(bounds.a.lower'*ones(1,N+1) < a < ...
        bounds.a.upper'*ones(1,N+1));
    opti.set_initial(a, guess.a');
    % Muscle activations at collocation points
    a_col = opti.variable(NMuscle,d*N);
    opti.subject_to(bounds.a.lower'*ones(1,d*N) < a_col < ...
        bounds.a.upper'*ones(1,d*N));
    opti.set_initial(a_col, guess.a_col');
    % Muscle-tendon forces at mesh points
    FTtilde = opti.variable(NMuscle,N+1);
    opti.subject_to(bounds.FTtilde.lower'*ones(1,N+1) < FTtilde < ...
        bounds.FTtilde.upper'*ones(1,N+1));
    opti.set_initial(FTtilde, guess.FTtilde');
    % Muscle-tendon forces at collocation points
    FTtilde_col = opti.variable(NMuscle,d*N);
    opti.subject_to(bounds.FTtilde.lower'*ones(1,d*N) < FTtilde_col < ...
        bounds.FTtilde.upper'*ones(1,d*N));
    opti.set_initial(FTtilde_col, guess.FTtilde_col');
    % Qs at mesh points
    Qs = opti.variable(nq.all,N+1);
    % We want to constraint the pelvis_tx position at the first mesh point,
    % and avoid redundant bounds
    lboundsQsk = bounds.QsQdots.lower(1:2:end)'*ones(1,N+1);
    lboundsQsk(jointi.pelvis.tx,1) = ...
        bounds.QsQdots_0.lower(2*jointi.pelvis.tx-1);
    uboundsQsk = bounds.QsQdots.upper(1:2:end)'*ones(1,N+1);
    uboundsQsk(jointi.pelvis.tx,1) = ...
        bounds.QsQdots_0.upper(2*jointi.pelvis.tx-1);
    opti.subject_to(lboundsQsk < Qs < uboundsQsk);
    opti.set_initial(Qs, guess.QsQdots(:,1:2:end)');
    % Qs at collocation points
    Qs_col = opti.variable(nq.all,d*N);
    opti.subject_to(bounds.QsQdots.lower(1:2:end)'*ones(1,d*N) < Qs_col < ...
        bounds.QsQdots.upper(1:2:end)'*ones(1,d*N));
    opti.set_initial(Qs_col, guess.QsQdots_col(:,1:2:end)');
    % Qdots at mesh points
    Qdots = opti.variable(nq.all,N+1);
    opti.subject_to(bounds.QsQdots.lower(2:2:end)'*ones(1,N+1) < Qdots < ...
        bounds.QsQdots.upper(2:2:end)'*ones(1,N+1));
    opti.set_initial(Qdots, guess.QsQdots(:,2:2:end)');
    % Qdots at collocation points
    Qdots_col = opti.variable(nq.all,d*N);
    opti.subject_to(bounds.QsQdots.lower(2:2:end)'*ones(1,d*N) < Qdots_col < ...
        bounds.QsQdots.upper(2:2:end)'*ones(1,d*N));
    opti.set_initial(Qdots_col, guess.QsQdots_col(:,2:2:end)');
    % Arm activations at mesh points
    a_a = opti.variable(nq.arms,N+1);
    opti.subject_to(bounds.a_a.lower'*ones(1,N+1) < a_a < ...
        bounds.a_a.upper'*ones(1,N+1));
    opti.set_initial(a_a, guess.a_a');
    % Arm activations at collocation points
    a_a_col = opti.variable(nq.arms,d*N);
    opti.subject_to(bounds.a_a.lower'*ones(1,d*N) < a_a_col < ...
        bounds.a_a.upper'*ones(1,d*N));
    opti.set_initial(a_a_col, guess.a_a_col');
    % Arm activations at mesh points
    a_mtp = opti.variable(nq.mtp,N+1);
    opti.subject_to(bounds.a_mtp.lower'*ones(1,N+1) < a_mtp < ...
        bounds.a_mtp.upper'*ones(1,N+1));
    opti.set_initial(a_mtp, guess.a_mtp');
    % Mtp activations at collocation points
    a_mtp_col = opti.variable(nq.mtp,d*N);
    opti.subject_to(bounds.a_mtp.lower'*ones(1,d*N) < a_mtp_col < ...
        bounds.a_mtp.upper'*ones(1,d*N));
    opti.set_initial(a_mtp_col, guess.a_mtp_col');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define controls
    % Time derivative of muscle activations (states) at mesh points
    vA = opti.variable(NMuscle, N);
    opti.subject_to(bounds.vA.lower'*ones(1,N) < vA < ...
        bounds.vA.upper'*ones(1,N));
    opti.set_initial(vA, guess.vA');
    % Arm excitations
    e_a = opti.variable(nq.arms, N);
    opti.subject_to(bounds.e_a.lower'*ones(1,N) < e_a < ...
        bounds.e_a.upper'*ones(1,N));
    opti.set_initial(e_a, guess.e_a');
    % Mtp excitations
    e_mtp = opti.variable(nq.mtp, N);
    opti.subject_to(bounds.e_mtp.lower'*ones(1,N) < e_mtp < ...
        bounds.e_mtp.upper'*ones(1,N));
    opti.set_initial(e_mtp, guess.e_mtp');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define "slack" controls
    % Time derivative of muscle-tendon forces (states) at collocation points
    dFTtilde_col = opti.variable(NMuscle, d*N);
    opti.subject_to(bounds.dFTtilde.lower'*ones(1,d*N) < dFTtilde_col < ...
        bounds.dFTtilde.upper'*ones(1,d*N));
    opti.set_initial(dFTtilde_col, guess.dFTtilde_col');
    % Time derivative of Qdots (states) at collocation points
    A_col = opti.variable(nq.all, d*N);
    opti.subject_to(bounds.Qdotdots.lower'*ones(1,d*N) < A_col < ...
        bounds.Qdotdots.upper'*ones(1,d*N));
    opti.set_initial(A_col, guess.Qdotdots_col');
    
    %% OCP: collocation equations
    % Define CasADi variables for static parameters
    tfk         = MX.sym('tfk');
    % Define CasADi variables for states
    ak          = MX.sym('ak',NMuscle);
    aj          = MX.sym('akmesh',NMuscle,d);
    akj         = [ak aj];
    FTtildek    = MX.sym('FTtildek',NMuscle);
    FTtildej    = MX.sym('FTtildej',NMuscle,d);
    FTtildekj   = [FTtildek FTtildej];
    Qsk         = MX.sym('Qsk',nq.all);
    Qsj         = MX.sym('Qsj',nq.all,d);
    Qskj        = [Qsk Qsj];
    Qdotsk      = MX.sym('Qdotsk',nq.all);
    Qdotsj      = MX.sym('Qdotsj',nq.all,d);
    Qdotskj     = [Qdotsk Qdotsj];
    a_ak        = MX.sym('a_ak',nq.arms);
    a_aj        = MX.sym('a_akmesh',nq.arms,d);
    a_akj       = [a_ak a_aj];
    a_mtpk      = MX.sym('a_mtpk',nq.mtp);
    a_mtpj      = MX.sym('a_mtpkmesh',nq.mtp,d);
    a_mtpkj     = [a_mtpk a_mtpj];
    
    % Define CasADi variables for controls
    vAk     = MX.sym('vAk',NMuscle);
    e_ak    = MX.sym('e_ak',nq.arms);
    e_mtpk  = MX.sym('e_mtpk',nq.mtp);
    
    % define the exoskeleton assistive torque
    Texok   = MX.sym('Texo',2); % joint moments for the exoskeleton
    
    % Define CasADi variables for "slack" controls
    dFTtildej   = MX.sym('dFTtildej',NMuscle,d);
    Aj          = MX.sym('Aj',nq.all,d);
    J           = 0; % Initialize cost function
    eq_constr   = {}; % Initialize equality constraint vector
    ineq_constr1 = {}; % Initialize inequality constraint vector 1
    ineq_constr2 = {}; % Initialize inequality constraint vector 2
    ineq_constr3 = {}; % Initialize inequality constraint vector 3
    ineq_constr4 = {}; % Initialize inequality constraint vector 4
    ineq_constr5 = {}; % Initialize inequality constraint vector 5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time step
    h = tfk/N;
    % Loop over collocation points
    for j=1:d
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Unscale variables
        Qskj_nsc = Qskj.*(scaling.QsQdots(1:2:end)'*ones(1,size(Qskj,2)/2));
        Qdotskj_nsc = Qdotskj.*(scaling.QsQdots(2:2:end)'* ...
            ones(1,size(Qdotskj,2)/2));
        FTtildekj_nsc = FTtildekj.*(scaling.FTtilde'*ones(1,size(FTtildekj,2)));
        dFTtildej_nsc = dFTtildej.*scaling.dFTtilde;
        Aj_nsc = Aj.*(scaling.Qdotdots'*ones(1,size(Aj,2)));
        vAk_nsc = vAk.*scaling.vA;
        
        QsQdotskj_nsc = MX(nq.all*2, d+1);
        QsQdotskj_nsc(1:2:end,:) = Qskj_nsc;
        QsQdotskj_nsc(2:2:end,:) = Qdotskj_nsc;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get muscle-tendon lengths, velocities, and moment arms
        % Left leg
        qinj_l          = Qskj_nsc(IndexLeft, j+1);
        qdotinj_l       = Qdotskj_nsc(IndexLeft, j+1);
        [lMTj_l,vMTj_l,MAj_l] =  f_lMT_vMT_dM(qinj_l,qdotinj_l);
        MAj.hip_flex.l   =  MAj_l(mai(1).mus.l',1);
        MAj.hip_add.l    =  MAj_l(mai(2).mus.l',2);
        MAj.hip_rot.l    =  MAj_l(mai(3).mus.l',3);
        MAj.knee.l       =  MAj_l(mai(4).mus.l',4);
        MAj.ankle.l      =  MAj_l(mai(5).mus.l',5);
        MAj.subt.l       =  MAj_l(mai(6).mus.l',6);
        % For the back muscles, we want left and right together: left
        % first, right second. In MuscleInfo, we first have the right
        % muscles (44:46) and then the left muscles (47:49). Since the back
        % muscles only depend on back dofs, we do not care if we extract
        % them "from the left or right leg" so here we just picked left.
        MAj.trunk_ext    =  MAj_l([47:49,mai(8).mus.l]',8);
        MAj.trunk_ben    =  MAj_l([47:49,mai(9).mus.l]',9);
        MAj.trunk_rot    =  MAj_l([47:49,mai(10).mus.l]',10);
        % Right leg
        qinj_r      = Qskj_nsc(IndexRight,j+1);
        qdotinj_r   = Qdotskj_nsc(IndexRight,j+1);
        [lMTj_r,vMTj_r,MAj_r] = f_lMT_vMT_dM(qinj_r,qdotinj_r);
        % Here we take the indices from left since the vector is 1:49
        MAj.hip_flex.r   =  MAj_r(mai(1).mus.l',1);
        MAj.hip_add.r    =  MAj_r(mai(2).mus.l',2);
        MAj.hip_rot.r    =  MAj_r(mai(3).mus.l',3);
        MAj.knee.r       =  MAj_r(mai(4).mus.l',4);
        MAj.ankle.r      =  MAj_r(mai(5).mus.l',5);
        MAj.subt.r       =  MAj_r(mai(6).mus.l',6);
        % Both legs
        % In MuscleInfo, we first have the right back muscles (44:46) and
        % then the left back muscles (47:49). Here we re-organize so that
        % we have first the left muscles and then the right muscles.
        lMTj_lr = [lMTj_l([1:43,47:49],1);lMTj_r(1:46,1)];
        vMTj_lr = [vMTj_l([1:43,47:49],1);vMTj_r(1:46,1)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get muscle-tendon forces and derive Hill-equilibrium
        [Hilldiffj,FTj,Fcej,Fpassj,Fisoj,vMmaxj,massMj] = ...
            f_forceEquilibrium_FtildeState_all_tendon(akj(:,j+1),...
            FTtildekj_nsc(:,j+1),dFTtildej_nsc(:,j),...
            lMTj_lr,vMTj_lr,tensions);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get metabolic energy rate if in the cost function
        % Get muscle fiber lengths
        [~,lMtildej] = f_FiberLength_TendonForce_tendon(...
            FTtildekj_nsc(:,j+1),lMTj_lr);
        % Get muscle fiber velocities
        [vMj,~] = f_FiberVelocity_TendonForce_tendon(...
            FTtildekj_nsc(:,j+1),dFTtildej_nsc(:,j),...
            lMTj_lr,vMTj_lr);
        % Get metabolic energy rate Bhargava et al. (2004)
        [e_totj,~,~,~,~,~] = fgetMetabolicEnergySmooth2004all(...
            akj(:,j+1),akj(:,j+1),lMtildej,...
            vMj,Fcej,Fpassj,massMj,pctsts,Fisoj,body_mass,10);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get passive joint torques
        
        Tau_passj_all = f_AllPassiveTorques(Qskj_nsc(:,j+1),Qdotskj_nsc(:,j+1));
        Tau_passj.hip.flex.l = Tau_passj_all(1);
        Tau_passj.hip.flex.r = Tau_passj_all(2);
        Tau_passj.hip.add.l = Tau_passj_all(3);
        Tau_passj.hip.add.r = Tau_passj_all(4);
        Tau_passj.hip.rot.l = Tau_passj_all(5);
        Tau_passj.hip.rot.r = Tau_passj_all(6);
        Tau_passj.knee.l = Tau_passj_all(7);
        Tau_passj.knee.r = Tau_passj_all(8);
        Tau_passj.ankle.l = Tau_passj_all(9);
        Tau_passj.ankle.r = Tau_passj_all(10);
        Tau_passj.subt.l = Tau_passj_all(11);
        Tau_passj.subt.r = Tau_passj_all(12);
        Tau_passj.mtp.all = Tau_passj_all(13:14);
        Tau_passj.trunk.ext = Tau_passj_all(15);
        Tau_passj.trunk.ben = Tau_passj_all(16);
        Tau_passj.trunk.rot = Tau_passj_all(17);
        Tau_passj.arm = Tau_passj_all(18:25);
        
        % Expression for the state derivatives at the collocation points
        Qsp_nsc      = Qskj_nsc*C(:,j+1);
        Qdotsp_nsc   = Qdotskj_nsc*C(:,j+1);
        FTtildep_nsc = FTtildekj_nsc*C(:,j+1);
        ap           = akj*C(:,j+1);
        a_ap         = a_akj*C(:,j+1);
        a_mtpp       = a_mtpkj*C(:,j+1);
        % Append collocation equations
        % Dynamic constraints are scaled using the same scale
        % factors as the ones used to scale the states
        % Activation dynamics (implicit formulation)
        eq_constr{end+1} = (h*vAk_nsc - ap)./scaling.a;
        % Contraction dynamics (implicit formulation)
        eq_constr{end+1} = (h*dFTtildej_nsc(:,j) - FTtildep_nsc)./...
            scaling.FTtilde';
        % Skeleton dynamics (implicit formulation)
        qdotj_nsc = Qdotskj_nsc(:,j+1); % velocity
        eq_constr{end+1} = (h*qdotj_nsc - Qsp_nsc)./scaling.QsQdots(1:2:end)';
        eq_constr{end+1} = (h*Aj_nsc(:,j) - Qdotsp_nsc)./...
            scaling.QsQdots(2:2:end)';
        % Arm activation dynamics (explicit formulation)
        da_adtj = f_ArmActivationDynamics(e_ak,a_akj(:,j+1)');
        eq_constr{end+1} = (h*da_adtj - a_ap)./scaling.a_a;
        % Mtp activation dynamics (explicit formulation)
        da_mtpdtj = f_MtpActivationDynamics(e_mtpk,a_mtpkj(:,j+1)');
        eq_constr{end+1} = (h*da_mtpdtj - a_mtpp);
        % Add contribution to the quadrature function
        J = J + 1*(...
            W.E*B(j+1)      *(f_J92exp(e_totj,exp_E))/body_mass*h + ...
            W.A*B(j+1)      *(f_J92(akj(:,j+1)'))*h + ...
            W.ArmE*B(j+1)   *(f_J8(e_ak))*h +...
            W.Mtp*B(j+1)    *(f_J2(e_mtpk))*h +...
            W.Ak*B(j+1)     *(f_J23(Aj(residuals_noarmsi,j)))*h + ...
            W.passMom*B(j+1)*(f_J25(Tau_passj_all))*h + ...
            W.u*B(j+1)      *(f_J92(vAk))*h + ...
            W.u*B(j+1)      *(f_J92(dFTtildej(:,j)))*h + ...
            W.u*B(j+1)      *(f_J8(Aj(armsi,j)))*h);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Call external function (run inverse dynamics)
        [Tj] = F([QsQdotskj_nsc(:,j+1);Aj_nsc(:,j)]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add path constraints
        % Null pelvis residuals
        eq_constr{end+1} = Tj(ground_pelvisi,1);
        % Muscle-driven joint torques for the lower limbs and the trunk
        % Hip flexion, left
        Ft_hip_flex_l   = FTj(mai(1).mus.l',1);
        T_hip_flex_l    = f_T27(MAj.hip_flex.l,Ft_hip_flex_l);
        eq_constr{end+1} = Tj(jointi.hip_flex.l,1)-(T_hip_flex_l + Tau_passj.hip.flex.l);
        % Hip flexion, right
        Ft_hip_flex_r   = FTj(mai(1).mus.r',1);
        T_hip_flex_r    = f_T27(MAj.hip_flex.r,Ft_hip_flex_r);
        eq_constr{end+1} = Tj(jointi.hip_flex.r,1)-(T_hip_flex_r + Tau_passj.hip.flex.r);
        % Hip adduction, left
        Ft_hip_add_l    = FTj(mai(2).mus.l',1);
        T_hip_add_l     = f_T27(MAj.hip_add.l,Ft_hip_add_l);
        eq_constr{end+1} = Tj(jointi.hip_add.l,1)-(T_hip_add_l + Tau_passj.hip.add.l);
        % Hip adduction, right
        Ft_hip_add_r    = FTj(mai(2).mus.r',1);
        T_hip_add_r     = f_T27(MAj.hip_add.r,Ft_hip_add_r);
        eq_constr{end+1} = Tj(jointi.hip_add.r,1)-(T_hip_add_r + Tau_passj.hip.add.r);
        % Hip rotation, left
        Ft_hip_rot_l    = FTj(mai(3).mus.l',1);
        T_hip_rot_l     = f_T27(MAj.hip_rot.l,Ft_hip_rot_l);
        eq_constr{end+1} = Tj(jointi.hip_rot.l,1)-(T_hip_rot_l + Tau_passj.hip.rot.l);
        % Hip rotation, right
        Ft_hip_rot_r    = FTj(mai(3).mus.r',1);
        T_hip_rot_r     = f_T27(MAj.hip_rot.r,Ft_hip_rot_r);
        eq_constr{end+1} = Tj(jointi.hip_rot.r,1)-(T_hip_rot_r + Tau_passj.hip.rot.r);
        % Knee, left
        Ft_knee_l       = FTj(mai(4).mus.l',1);
        T_knee_l        = f_T13(MAj.knee.l,Ft_knee_l);
        eq_constr{end+1} = Tj(jointi.knee.l,1)-(T_knee_l + Tau_passj.knee.l);
        % Knee, right
        Ft_knee_r       = FTj(mai(4).mus.r',1);
        T_knee_r        = f_T13(MAj.knee.r,Ft_knee_r);
        eq_constr{end+1} = Tj(jointi.knee.r,1)-(T_knee_r + Tau_passj.knee.r);
        % Ankle, left
        Ft_ankle_l      = FTj(mai(5).mus.l',1);
        T_ankle_l       = f_T12(MAj.ankle.l,Ft_ankle_l);
        eq_constr{end+1} = Tj(jointi.ankle.l,1)-(T_ankle_l + Tau_passj.ankle.l  + Texok(1));
        % Ankle, right
        Ft_ankle_r      = FTj(mai(5).mus.r',1);
        T_ankle_r       = f_T12(MAj.ankle.r,Ft_ankle_r);
        eq_constr{end+1} = Tj(jointi.ankle.r,1)-(T_ankle_r + Tau_passj.ankle.r + Texok(2));
        % Subtalar, left
        Ft_subt_l       = FTj(mai(6).mus.l',1);
        T_subt_l        = f_T12(MAj.subt.l,Ft_subt_l);
        eq_constr{end+1} = Tj(jointi.subt.l,1)-(T_subt_l +  Tau_passj.subt.l);
        % Subtalar, right
        Ft_subt_r       = FTj(mai(6).mus.r',1);
        T_subt_r        = f_T12(MAj.subt.r,Ft_subt_r);
        eq_constr{end+1} = Tj(jointi.subt.r,1)-(T_subt_r + Tau_passj.subt.r );
        % Lumbar extension
        Ft_trunk_ext    = FTj([mai(8).mus.l,mai(8).mus.r]',1);
        T_trunk_ext     = f_T6(MAj.trunk_ext,Ft_trunk_ext);
        eq_constr{end+1} = Tj(jointi.trunk.ext,1)-(T_trunk_ext + Tau_passj.trunk.ext);
        % Lumbar bending
        Ft_trunk_ben    = FTj([mai(9).mus.l,mai(9).mus.r]',1);
        T_trunk_ben     = f_T6(MAj.trunk_ben,Ft_trunk_ben);
        eq_constr{end+1} = Tj(jointi.trunk.ben,1)-(T_trunk_ben + Tau_passj.trunk.ben);
        % Lumbar rotation
        Ft_trunk_rot    = FTj([mai(10).mus.l,mai(10).mus.r]',1);
        T_trunk_rot     = f_T6(MAj.trunk_rot,Ft_trunk_rot);
        eq_constr{end+1} = Tj(jointi.trunk.rot,1)-(T_trunk_rot + Tau_passj.trunk.rot);
        % Torque-driven joint torques for the arms
        % Arms
        eq_constr{end+1} = Tj(armsi,1)/scaling.ArmTau - (a_akj(:,j+1) + ...
            (Tau_passj.arm)/scaling.ArmTau);
        % Mtp
        eq_constr{end+1} = Tj(mtpi,1)/scaling.MtpTau - (a_mtpkj(:,j+1) + ...
            (Tau_passj.mtp.all)/scaling.MtpTau);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Activation dynamics (implicit formulation)
        act1 = vAk_nsc + akj(:,j+1)./(ones(size(akj(:,j+1),1),1)*tdeact);
        act2 = vAk_nsc + akj(:,j+1)./(ones(size(akj(:,j+1),1),1)*tact);
        ineq_constr1{end+1} = act1;
        ineq_constr2{end+1} = act2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Contraction dynamics (implicit formulation)
        eq_constr{end+1} = Hilldiffj;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constraints to prevent parts of the skeleton to penetrate each
        % other.
        % Origins calcaneus (transv plane) at minimum 9 cm from each other.
        ineq_constr3{end+1} = f_Jnn2(Tj(calcOr.r,1) - Tj(calcOr.l,1));
        % Constraint to prevent the arms to penetrate the skeleton
        % Origins femurs and ipsilateral hands (transv plane) at minimum
        % 18 cm from each other.
        ineq_constr4{end+1} = f_Jnn2(Tj(femurOr.r,1) - Tj(handOr.r,1));
        ineq_constr4{end+1} = f_Jnn2(Tj(femurOr.l,1) - Tj(handOr.l,1));
        % Origins tibia (transv plane) at minimum 11 cm from each other.
        ineq_constr5{end+1} = f_Jnn2(Tj(tibiaOr.r,1) - Tj(tibiaOr.l,1));
    end % End loop over collocation points
    eq_constr = vertcat(eq_constr{:});
    ineq_constr1 = vertcat(ineq_constr1{:});
    ineq_constr2 = vertcat(ineq_constr2{:});
    ineq_constr3 = vertcat(ineq_constr3{:});
    ineq_constr4 = vertcat(ineq_constr4{:});
    ineq_constr5 = vertcat(ineq_constr5{:});
    f_coll = Function('f_coll',{tfk,ak,aj,FTtildek,FTtildej,Qsk,Qsj,Qdotsk,...
        Qdotsj,a_ak,a_aj,a_mtpk,a_mtpj,vAk,e_ak,e_mtpk,dFTtildej,Aj,Texok},...
        {eq_constr,ineq_constr1,ineq_constr2,ineq_constr3,ineq_constr4,...
        ineq_constr5,J});
    f_coll_map = f_coll.map(N,parallelMode,NThreads);
    [coll_eq_constr, coll_ineq_constr1, coll_ineq_constr2, coll_ineq_constr3,...
        coll_ineq_constr4, coll_ineq_constr5, Jall] = f_coll_map(tf,...
        a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, Qs(:,1:end-1), ...
        Qs_col, Qdots(:,1:end-1), Qdots_col, a_a(:,1:end-1), a_a_col, ...
        a_mtp(:,1:end-1), a_mtp_col, vA, e_a, e_mtp, dFTtilde_col, A_col,ExoVect);
    opti.subject_to(coll_eq_constr == 0);
    opti.subject_to(coll_ineq_constr1(:) >= 0);
    opti.subject_to(coll_ineq_constr2(:) <= 1/tact);
    opti.subject_to(0.0081 < coll_ineq_constr3(:) < 4);
    opti.subject_to(0.0324 < coll_ineq_constr4(:) < 4);
    opti.subject_to(0.0121 < coll_ineq_constr5(:) < 4);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop over mesh points
    for k=1:N
        % Variables within current mesh interval
        % States
        akj = [a(:,k), a_col(:,(k-1)*d+1:k*d)];
        FTtildekj = [FTtilde(:,k), FTtilde_col(:,(k-1)*d+1:k*d)];
        Qskj = [Qs(:,k), Qs_col(:,(k-1)*d+1:k*d)];
        Qdotskj = [Qdots(:,k), Qdots_col(:,(k-1)*d+1:k*d)];
        a_akj = [a_a(:,k), a_a_col(:,(k-1)*d+1:k*d)];
        a_mtpkj = [a_mtp(:,k), a_mtp_col(:,(k-1)*d+1:k*d)];
        % Add equality constraints (next interval starts with end values of
        % states from previous interval)
        opti.subject_to(a(:,k+1) == akj*D);
        opti.subject_to(FTtilde(:,k+1) == FTtildekj*D); % scaled
        opti.subject_to(Qs(:,k+1) == Qskj*D); % scaled
        opti.subject_to(Qdots(:,k+1) == Qdotskj*D); % scaled
        opti.subject_to(a_a(:,k+1) == a_akj*D);
        opti.subject_to(a_mtp(:,k+1) == a_mtpkj*D);
    end % End loop over mesh points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Additional path constraints
    % Periodicity of the states
    % Qs and Qdots
    
    opti.subject_to(Qs(QsInvA,end) - Qs(QsInvB,1) == 0);
    opti.subject_to(Qdots(QdotsInvA,end) - Qdots(QdotsInvB,1) == 0);
    opti.subject_to(Qs(orderQsOpp,end) + Qs(orderQsOpp,1) == 0);
    opti.subject_to(Qdots(orderQsOpp,end) + Qdots(orderQsOpp,1) == 0);
    % Muscle activations
    orderMusInv = [NMuscle/2+1:NMuscle,1:NMuscle/2];
    opti.subject_to(a(:,end) - a(orderMusInv,1) == 0);
    % Muscle-tendon forces
    opti.subject_to(FTtilde(:,end) - FTtilde(orderMusInv,1) == 0);
    % Arm activations
    orderArmInv = [jointi.sh_flex.r:jointi.sh_rot.r,...
        jointi.sh_flex.l:jointi.sh_rot.l,...
        jointi.elb.r:jointi.elb.r,...
        jointi.elb.l:jointi.elb.l]-jointi.sh_flex.l+1;
    opti.subject_to(a_a(:,end) - a_a(orderArmInv,1) == 0);
    % Mtp activations
    orderMtpInv = [jointi.mtp.r,jointi.mtp.l]-jointi.mtp.l+1;
    opti.subject_to(a_mtp(:,end) - a_mtp(orderMtpInv,1) == 0);
    % Average speed
    % Provide expression for the distance traveled
    Qs_nsc = Qs.*(scaling.QsQdots(1:2:end)'*ones(1,N+1));
    dist_trav_tot = Qs_nsc(jointi.pelvis.tx,end) -  Qs_nsc(jointi.pelvis.tx,1);
    vel_aver_tot = dist_trav_tot/tf;
    opti.subject_to(vel_aver_tot - v_tgt == 0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scale cost function
    Jall_sc = sum(Jall)/dist_trav_tot;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create NLP solver
    opti.minimize(Jall_sc);
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy           = 'adaptive';
    options.ipopt.max_iter              = 10000;
    options.ipopt.linear_solver         = linear_solver;
    options.ipopt.tol                   = 1*10^(-tol_ipopt);
    opti.solver('ipopt', options);
    % Create and save diary
    Outname = fullfile(OutFolder,[S.savename '_log.txt']);
    diary(Outname);
    % Data-informed (full solution at closest speed) initial guess
    if IGm == 4
        disp('Not supported')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve problem
    % Opti does not use bounds on variables but constraints. This function
    % adjusts for that.
    %     opti.solve();
    [w_opt,stats] = solve_NLPSOL(opti,options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    diary off
    % Extract results
    % Create setup
    setup.tolerance.ipopt = tol_ipopt;
    setup.bounds = bounds;
    setup.scaling = scaling;
    setup.guess = guess;
    
    Outname = fullfile(OutFolder,[S.savename '.mat']);
    Sopt = S;
    save(Outname,'w_opt','stats','setup','Sopt','ExoControl');    
end


end

