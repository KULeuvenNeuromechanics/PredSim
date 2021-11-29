function [] = f_PredSim_Gait92(S)

%% Adding the casadi path seems to be needed to run processes in batch
AddCasadiPaths();

%% Default settings
S = GetDefaultSettings(S);

%% User inputs (typical settings structure)
% settings for optimization
N           = S.N;          % number of mesh intervals
W           = S.W;          % weights optimization

%% Load external functions
import casadi.*
% The external function performs inverse dynamics through the
% OpenSim/Simbody C++ API. This external function is compiled as a dll from
% which we create a Function instance using CasADi in MATLAB. More details
% about the external function can be found in the documentation.
pathmain = pwd;
[filepath,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(filepath);
addpath(genpath(pathRepo));
% Loading external functions.
setup.derivatives =  'AD'; % Algorithmic differentiation
pathExternalFunctions = [pathRepo,'/ExternalFunctions'];
cd(pathExternalFunctions)
F  = external('F',S.ExternalFunc);
cd(pathmain);

%% Indices external function
% Indices of the elements in the external functions
% External function: F
% First, joint torques.
jointi = getJointi();

% Vectors of indices for later use
residualsi          = jointi.pelvis.tilt:jointi.elb.r; % all
ground_pelvisi      = jointi.pelvis.tilt:jointi.pelvis.tz; % ground-pelvis
trunki              = jointi.trunk.ext:jointi.trunk.rot; % trunk
armsi               = jointi.sh_flex.l:jointi.elb.r; % arms
mtpi                = jointi.mtp.l:jointi.mtp.r; % mtps
residuals_noarmsi   = jointi.pelvis.tilt:jointi.trunk.rot; % all but arms
% Number of degrees of freedom for later use
nq.all      = length(residualsi); % all
nq.abs      = length(ground_pelvisi); % ground-pelvis
nq.trunk    = length(trunki); % trunk
nq.arms     = length(armsi); % arms
nq.mtp      = length(mtpi); % arms
nq.leg      = 10; % #joints needed for polynomials
% Second, origins bodies.
% Calcaneus
calcOr.r    = 32:33;
calcOr.l    = 34:35;
calcOr.all  = [calcOr.r,calcOr.l];
% Femurs
femurOr.r   = 36:37;
femurOr.l   = 38:39;
femurOr.all = [femurOr.r,femurOr.l];
% Hands
handOr.r    = 40:41;
handOr.l    = 42:43;
handOr.all  = [handOr.r,handOr.l];
% Tibias
tibiaOr.r   = 44:45;
tibiaOr.l   = 46:47;
tibiaOr.all = [tibiaOr.r,tibiaOr.l];
% toe joints
toesOr.r   = 48:49;
toesOr.l   = 50:51;
toesOr.all = [toesOr.r,toesOr.l];
% velocity calcaneus 
calcOr_dot.r = 52:54;
calcOr_dot.l = 55:57;
calcOr_dot.all = [calcOr_dot.r calcOr_dot.l];
% velocity toes
toesOr_dot.r   = 58:60;
toesOr_dot.l   = 61:63;
toesOr_dot.all = [toesOr_dot.r,toesOr_dot.l];

calcOry.r    = 64;
calcOry.l    = 65;
calcOry.all  = [calcOry.r,calcOry.l];

toesOry.r    = 66;
toesOry.l    = 67;
toesOry.all  = [toesOry.r, toesOry.l];



%% Collocation scheme
% We use a pseudospectral direct collocation method, i.e. we use Lagrange
% polynomials to approximate the state derivatives at the collocation
% points in each mesh interval. We use d=3 collocation points per mesh
% interval and Radau collocation points.
pathCollocationScheme = [pathRepo,'/CollocationScheme'];
addpath(genpath(pathCollocationScheme));
d = 3; % degree of interpolating polynomial
method = 'radau'; % collocation method
[~,C,D,B] = CollocationScheme(d,method);

%% Muscle information
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
% Total number of muscles
NMuscle = length(muscleNames(1:end-3))*2;
pathpolynomial = fullfile(pathRepo,'Polynomials',S.PolyFolder); % default location
tl = load([pathpolynomial,'/muscle_spanning_joint_INFO.mat']);
[~,mai] = MomentArmIndices(muscleNames(1:end-3),tl.muscle_spanning_joint_INFO);

% Parameters for activation dynamics
tact = 0.015; % Activation time constant
tdeact = 0.06; % Deactivation time constant

%% Metabolic energy model parameters
% We extract the specific tensions and slow twitch rations.
pathMetabolicEnergy = [pathRepo,'/MetabolicEnergy'];
addpath(genpath(pathMetabolicEnergy));
tension = getSpecificTensions(muscleNames(1:end-3));
tensions = [tension;tension];
pctst = getSlowTwitchRatios(muscleNames(1:end-3));
pctsts = [pctst;pctst];

%% CasADi functions
% We create several CasADi functions for later use
pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
PathDefaultFunc = fullfile(pathCasADiFunctions,S.CasadiFunc_Folders);
f_ArmActivationDynamics = Function.load(fullfile(PathDefaultFunc,'f_ArmActivationDynamics'));
f_FiberLength_TendonForce_tendon = Function.load(fullfile(PathDefaultFunc,'f_FiberLength_TendonForce_tendon'));
f_FiberVelocity_TendonForce_tendon = Function.load(fullfile(PathDefaultFunc,'f_FiberVelocity_TendonForce_tendon'));
f_forceEquilibrium_FtildeState_all_tendon = Function.load(fullfile(PathDefaultFunc,'f_forceEquilibrium_FtildeState_all_tendon'));
f_J2    = Function.load(fullfile(PathDefaultFunc,'f_J2'));
f_J23   = Function.load(fullfile(PathDefaultFunc,'f_J23'));
f_J8    = Function.load(fullfile(PathDefaultFunc,'f_J8'));
f_J92   = Function.load(fullfile(PathDefaultFunc,'f_J92'));
f_J92exp = Function.load(fullfile(PathDefaultFunc,'f_J92exp'));
f_Jnn2  = Function.load(fullfile(PathDefaultFunc,'f_Jnn2'));
f_lMT_vMT_dM = Function.load(fullfile(PathDefaultFunc,'f_lMT_vMT_dM'));
f_MtpActivationDynamics = Function.load(fullfile(PathDefaultFunc,'f_MtpActivationDynamics'));
f_T12 = Function.load(fullfile(PathDefaultFunc,'f_T12'));
f_T13 = Function.load(fullfile(PathDefaultFunc,'f_T13'));
f_T27 = Function.load(fullfile(PathDefaultFunc,'f_T27'));
f_T6 = Function.load(fullfile(PathDefaultFunc,'f_T6'));
f_AllPassiveTorques = Function.load(fullfile(PathDefaultFunc,'f_AllPassiveTorques'));
fgetMetabolicEnergySmooth2004all = Function.load(fullfile(PathDefaultFunc,'fgetMetabolicEnergySmooth2004all'));

PathEnergyEq = fullfile(pathCasADiFunctions,'EnergyModels');
fgetMetabolicEnergy_MargariaSmooth  = Function.load(fullfile(PathEnergyEq, 'fgetMetabolicEnergy_MargariaSmooth'));


% file with mass of muscles
MassFile = fullfile(PathDefaultFunc,'MassM.mat');
if exist(MassFile,'file')
    MuscleMass = load(MassFile);
else
    MassFile = fullfile(pathCasADiFunctions,'MassM.mat'); % default muscle mass
    MuscleMass =load(MassFile);
end



%% Get bounds and initial guess

% Kinematics file for bounds -- input arguments
IKfile_bounds = fullfile(pathRepo, S.IKfile_Bounds);

% We extract experimental data to set bounds and initial guesses if needed
joints = {'pelvis_tilt','pelvis_list','pelvis_rotation','pelvis_tx',...
    'pelvis_ty','pelvis_tz','hip_flexion_l','hip_adduction_l',...
    'hip_rotation_l','hip_flexion_r','hip_adduction_r','hip_rotation_r',...
    'knee_angle_l','knee_angle_r','ankle_angle_l','ankle_angle_r',...
    'subtalar_angle_l','subtalar_angle_r','mtp_angle_l','mtp_angle_r',...
    'lumbar_extension','lumbar_bending','lumbar_rotation','arm_flex_l',...
    'arm_add_l','arm_rot_l','arm_flex_r','arm_add_r','arm_rot_r',...
    'elbow_flex_l','elbow_flex_r'};
Qs_walk          = getIK(IKfile_bounds,joints);

if S.Bounds_Running 
    [bounds,scaling] = getBounds_RunningMaartenB(Qs_walk,NMuscle,nq,jointi,S.v_tgt);
else
    [bounds,scaling] = getBounds_all_mtp(Qs_walk,NMuscle,nq,jointi,S.v_tgt);
end

% adapt bounds based on user input
bounds = AdaptBounds(bounds,S,mai);

% The initial guess depends on the settings
pathIG = [pathRepo,'/IG'];
addpath(genpath(pathIG));
if S.IGsel == 1 % Quasi-random initial guess
    guess = getGuess_QR_opti_int(N,nq,NMuscle,scaling,S.v_tgt,jointi,d,S.IG_PelvisY);
elseif S.IGsel == 2 % Data-informed initial guess
    if S.IGmodeID  < 3 % Data from average walking motion
        IKfile_guess    = fullfile(pathRepo, S.IKfile_guess);
        Qs_guess        = getIK(IKfile_guess,joints);
        time_IC         = [Qs_guess.time(1),Qs_guess.time(end)];
        guess = getGuess_DI_opti_int_mtp(Qs_guess,nq,N,time_IC,NMuscle,jointi,...
            scaling,S.v_tgt,d);
    elseif S.IGmodeID == 3 || S.IGmodeID == 4 % Data from selected motion
        % Extract joint positions from existing motion (previous results)
        if S.IGmodeID == 3
            GuessFolder = fullfile(pathRepo,'Results',S.ResultsF_ig);
        elseif S.IGmodeID ==4
            GuessFolder = fullfile(pathRepo,'IG','data');
        end
        pathIK      = fullfile(GuessFolder,[S.savename_ig '.mot']);
        Qs_ig       = getIK(pathIK,joints);
        % When saving the results, we save a 2 full gait cycle (4*N) so here we
        % only select 1:N to have half a gait cycle
        nfr = length(Qs_ig.allfilt(:,1));
        frSel = round(nfr./4);
        Qs_ig_sel.allfilt   = Qs_ig.allfilt(1:frSel,:);
        Qs_ig_sel.time      = Qs_ig.time(1:frSel,:);
        Qs_ig_sel.colheaders = Qs_ig.colheaders;
        time_IC = [Qs_ig_sel.time(1),Qs_ig_sel.time(end)];
        guess = getGuess_DI_opti_int_mtp(Qs_ig_sel,nq,N,time_IC,NMuscle,jointi,scaling,S.v_tgt,d);
    end
end

% adapt guess so that it fits within the bounds
guess = AdaptGuess_UserInput(guess,bounds,S);

%% exoskeleton torques
% function to get exoskeleton torques at mesh points
if ~S.ExoDesignBool
    ExoVect = GetExoTorques(S,pathRepo,N);
end

%% Index helpers
% get help indexes for left and right leg and for symmetry constraint
[IndexLeft,IndexRight,QsInvA,QsInvB,QdotsInvA,...
    QdotsInvB,orderQsOpp] = GetIndexHelper(S,jointi);

%% OCP create variables and bounds
% using opti
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design optimal exoskeleton support
if S.OptTexo_Ankle.Bool
    ExoVect = opti.variable(2, N);
    opti.subject_to(S.OptTexo_Ankle.Tbound(1)< ExoVect < S.OptTexo_Ankle.Tbound(2));
    opti.set_initial(ExoVect,zeros(2,N));
elseif S.OptTexo_AnkleKneeHip.Bool
    ExoVect = opti.variable(6, N);
    opti.subject_to(S.OptTexo_AnkleKneeHip.Tbound_Ankle(1) < ...
        ExoVect(1:2,:) < S.OptTexo_AnkleKneeHip.Tbound_Ankle(2));
    opti.subject_to(S.OptTexo_AnkleKneeHip.Tbound_Knee(1) < ...
        ExoVect(3:4,:) < S.OptTexo_AnkleKneeHip.Tbound_Knee(2));
    opti.subject_to(S.OptTexo_AnkleKneeHip.Tbound_Hip(1) < ...
        ExoVect(5:6,:) < S.OptTexo_AnkleKneeHip.Tbound_Hip(2));
end

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
nExoDofs = length(ExoVect(:,1));
Texok   = MX.sym('Texo',nExoDofs,1); % joint moments for the exoskeleton

% Define CasADi variables for "slack" controls
dFTtildej   = MX.sym('dFTtildej',NMuscle,d);
Aj          = MX.sym('Aj',nq.all,d);
J           = 0; % Initialize cost function
eq_constr   = {}; % Initialize equality constraint vector
ineq_constr = MX(1000,1); % Initialise inequality constraint vector
ineq_constr_index = nan(1000); % keeps track of the index of the inequality constraint
ctIneq      = 1;
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
    [Hilldiffj,FTj,Fcej,Fpassj,Fisoj] = ...
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
    if strcmp(S.EModel,'Bhargava2014')
        % Get metabolic energy rate Bhargava et al. (2004)
        [e_totj,~,~,~,~,~] = fgetMetabolicEnergySmooth2004all(...
            akj(:,j+1),akj(:,j+1),lMtildej,vMj,Fcej,Fpassj,...
            MuscleMass.MassM',pctsts,Fisoj,S.mass,10);
    elseif  strcmp(S.EModel,'Marg1968')
        e_totj = fgetMetabolicEnergy_MargariaSmooth(Fcej,vMj);
    else
        error('No energy model selected');
    end
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
    Tau_passj_J = Tau_passj_all([1:12 15:end]);
    
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
        W.E*B(j+1)      *(f_J92exp(e_totj,W.exp_E))/S.mass*h + ...
        W.A*B(j+1)      *(f_J92(akj(:,j+1)'))*h + ...
        W.ArmE*B(j+1)   *(f_J8(e_ak))*h +...
        W.Mtp*B(j+1)    *(f_J2(e_mtpk))*h +...
        W.Ak*B(j+1)     *(f_J23(Aj(residuals_noarmsi,j)))*h + ...
        W.passMom*B(j+1)*(f_J23(Tau_passj_J))*h + ...
        W.u*B(j+1)      *(f_J92(vAk))*h + ...
        W.u*B(j+1)      *(f_J92(dFTtildej(:,j)))*h + ...
        W.u*B(j+1)      *(f_J8(Aj(armsi,j)))*h);

    if W.EccF ~= 0
        % select eccentric muscle force
        b = 100;
        ivMPos = (tanh(b*(vMj-0.02))*0.5+0.5);      % negative power = 0, positive power = positive power
        Fecc = FTj.*ivMPos;
        % add Fecc squared to the objective function
        J = J + 1*(...
            W.EccF*B(j+1)    *(f_J92exp(Fecc,2))*h );
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Call external function (run inverse dynamics)
    if F.nnz_in == nq.all*3
        % no exoskeleton torque as input in passive simulations
        [Tj] = F([QsQdotskj_nsc(:,j+1);Aj_nsc(:,j)]);    % left and right leg exoskeleton torques as inputs as well.
    elseif F.nnz_in > nq.all*3
        % exoskeleton torques as input in active simulations
        [Tj] = F([QsQdotskj_nsc(:,j+1);Aj_nsc(:,j); -Texok]);    % left and right leg exoskeleton torques as inputs as well.
    end
    
    % note that this has to be -Texo, since a positive torque around
    % the z-axis of the tibia results in a plantarflexion moment. (and
    % plantarflexion is defined as negative in opensim and in the
    % exo torque vector.
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
    eq_constr{end+1} = Tj(jointi.ankle.l,1)-(T_ankle_l + Tau_passj.ankle.l);
    % Ankle, right
    Ft_ankle_r      = FTj(mai(5).mus.r',1);
    T_ankle_r       = f_T12(MAj.ankle.r,Ft_ankle_r);
    eq_constr{end+1} = Tj(jointi.ankle.r,1)-(T_ankle_r + Tau_passj.ankle.r);
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
    ineq_constr(ctIneq:ctIneq+length(act1)-1) = act1;
    ineq_constr_index(ctIneq:ctIneq+length(act1)-1) = 1;
    ctIneq = ctIneq + length(act1);
    ineq_constr(ctIneq:ctIneq+length(act2)-1) = act2;
    ineq_constr_index(ctIneq:ctIneq+length(act2)-1) = 2;
    ctIneq = ctIneq + length(act2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contraction dynamics (implicit formulation)
    eq_constr{end+1} = Hilldiffj;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constraints to prevent parts of the skeleton to penetrate each
    % other.
    % Origins calcaneus (transv plane) at minimum 9 cm from each other.
    Qconstr = f_Jnn2(Tj(calcOr.r,1) - Tj(calcOr.l,1));
    ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) = Qconstr;
    ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 3;
    ctIneq = ctIneq + length(Qconstr);
    % Constraint to prevent the arms to penetrate the skeleton
    % Origins femurs and ipsilateral hands (transv plane) at minimum
    % 18 cm from each other.
    Qconstr = f_Jnn2(Tj(femurOr.r,1) - Tj(handOr.r,1));
    ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) = Qconstr;
    ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 4;
    ctIneq = ctIneq + length(Qconstr);
    Qconstr = f_Jnn2(Tj(femurOr.l,1) - Tj(handOr.l,1));
    ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) = Qconstr;
    ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 4;
    ctIneq = ctIneq + length(Qconstr);
    % Origins tibia (transv plane) at minimum 11 cm from each other.
    Qconstr = f_Jnn2(Tj(tibiaOr.r,1) - Tj(tibiaOr.l,1));
    ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) = Qconstr;
    ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 5;
    ctIneq = ctIneq + length(Qconstr);
    % Origins toes (transv plane) at minimum 10 cm from each other.
    Qconstr = f_Jnn2(Tj(toesOr.r,1) - Tj(toesOr.l,1));
    ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) = Qconstr;
    ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 6;
    ctIneq = ctIneq + length(Qconstr);
    % Minimum foot height during swing
    if isfield(S.Constr,'ImposeFootHeight') && S.Constr.ImposeFootHeight
        Vx = Tj([toesOr_dot.r(1) toesOr_dot.l(1) calcOr_dot.r(1) calcOr_dot.l(1)]);
        Ry = Tj([toesOry.r toesOry.l calcOry.r calcOry.l]);
        vx_tr = S.Constr.FootHeight.Vx_Treshold;
        tr_ry = S.Constr.FootHeight.Ry_Treshold;
        b = S.Constr.FootHeight.b;
        vx_bool = 0.5.*tanh(b*(Vx-vx_tr))+0.5;
        ry_ad = -Ry+tr_ry;
        Qconstr =  ry_ad.*vx_bool;
        ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) =Qconstr;
        ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 101;
        ctIneq = ctIneq + length(Qconstr);
    end
    % inequality constraint on exoskeleton moments
    if S.OptTexo_Ankle.Bool
        Qconstr = Texok(1).*Qskj_nsc(jointi.ankle.l,j+1);
        ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) = Qconstr;
        ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 7;
        ctIneq = ctIneq + length(Qconstr);
        Qconstr = Texok(2).*Qskj_nsc(jointi.ankle.r,j+1);
        ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) = Qconstr;
        ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 7;
        ctIneq = ctIneq + length(Qconstr);
    elseif S.OptTexo_AnkleKneeHip.Bool
        iVect = [jointi.ankle.l jointi.ankle.r jointi.knee.l jointi.knee.r , ...
            jointi.hip_flex.l jointi.hip_flex.r];
        for iv = 1:6
            Qconstr = Texok(iv).*Qskj_nsc(iVect(iv),j+1);
            ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) = Qconstr;
            ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 7+iv-1;
            ctIneq = ctIneq + length(Qconstr);
        end        
    end
end % End loop over collocation points
eq_constrV = vertcat(eq_constr{:});

% Casadi function to get constraints and objective
f_coll = Function('f_coll',{tfk,ak,aj,FTtildek,FTtildej,Qsk,Qsj,Qdotsk,...
    Qdotsj,a_ak,a_aj,a_mtpk,a_mtpj,vAk,e_ak,e_mtpk,dFTtildej,Aj,Texok},...
    {eq_constrV,ineq_constr,J});
% assign NLP problem to multiple cores
f_coll_map = f_coll.map(N,S.parallelMode,S.NThreads);
[coll_eq_constr, coll_ineq_constr, Jall] = f_coll_map(tf,...
    a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, Qs(:,1:end-1), ...
    Qs_col, Qdots(:,1:end-1), Qdots_col, a_a(:,1:end-1), a_a_col, ...
    a_mtp(:,1:end-1), a_mtp_col, vA, e_a, e_mtp, dFTtilde_col, A_col,ExoVect);
% equality constrains
opti.subject_to(coll_eq_constr == 0);

% inequality constraints (logical indexing not possible in MX arrays)
cSel = coll_ineq_constr(find(ineq_constr_index == 1),:); % activation dyanmics
opti.subject_to(cSel(:)  >= 0);
cSel = coll_ineq_constr(find(ineq_constr_index == 2),:); % deactivation dyanmics
opti.subject_to(cSel(:)  <= 1/tact);
cSel = coll_ineq_constr(find(ineq_constr_index == 3),:); % origin calcaneus
opti.subject_to(S.Constr.calcn.^2 < cSel(:) < 4);
cSel = coll_ineq_constr(find(ineq_constr_index == 4),:); % arms
opti.subject_to(0.0324 < cSel(:) < 4);
cSel = coll_ineq_constr(find(ineq_constr_index == 5),:); % origin tibia minimum x cm away from each other
opti.subject_to(S.Constr.tibia.^2 < cSel(:) < 4);
cSel = coll_ineq_constr(find(ineq_constr_index == 6),:); % origins toes minimum x cm away from each other
opti.subject_to(S.Constr.toes.^2 < cSel(:) < 4);
if S.OptTexo_Ankle.Bool         % bound on exoskeleton power
    cSel = coll_ineq_constr(find(ineq_constr_index == 7),:); % exoskeleton power
    opti.subject_to(S.OptTexo_Ankle.Pbound(1) < cSel(:) < S.OptTexo_Ankle.Pbound(2));
elseif S.OptTexo_AnkleKneeHip.Bool
    lb_ExoPower = [S.OptTexo_AnkleKneeHip.Pbound_Ankle(1)*ones(1,2) ,...
        S.OptTexo_AnkleKneeHip.Pbound_Knee(1)*ones(1,2), ...
        S.OptTexo_AnkleKneeHip.Pbound_Hip(1)*ones(1,2)];
    ub_ExoPower = [S.OptTexo_AnkleKneeHip.Pbound_Ankle(2)*ones(1,2) ,...
        S.OptTexo_AnkleKneeHip.Pbound_Knee(2)*ones(1,2), ...
        S.OptTexo_AnkleKneeHip.Pbound_Hip(2)*ones(1,2)];
    for iv = 1:6
        cSel = coll_ineq_constr(find(ineq_constr_index == 7+iv-1),:); % exoskeleton power
        opti.subject_to(lb_ExoPower(iv) < cSel(:) < ub_ExoPower(iv));
    end
end
% inequality constraint minimum foot height during swing
if isfield(S.Constr,'ImposeFootHeight') && S.Constr.ImposeFootHeight
    cSel = coll_ineq_constr(find(ineq_constr_index == 101),:); % origins toes minimum x cm away from each other
    opti.subject_to(cSel(:) < 0.01);
end

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
if S.Symmetric
    % Periodicity of the states (or rather LR symmetry -half gait cycle)
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
elseif S.Periodic
    opti.subject_to(Qs(:,end) - Qs(:,1) == 0);
    opti.subject_to(Qdots(:,end) - Qdots(:,1) == 0);
    opti.subject_to(Qs(:,end) + Qs(:,1) == 0);
    opti.subject_to(Qdots(:,end) + Qdots(:,1) == 0);
    % Muscle activations
    opti.subject_to(a(:,end) - a(:,1) == 0);
    % Muscle-tendon forces
    opti.subject_to(FTtilde(:,end) - FTtilde(:,1) == 0);
    % Arm activations
    opti.subject_to(a_a(:,end) - a_a(:,1) == 0);
    % Mtp activations
    opti.subject_to(a_mtp(:,end) - a_mtp(:,1) == 0);
end
% Average speed
% Provide expression for the distance traveled
Qs_nsc = Qs.*(scaling.QsQdots(1:2:end)'*ones(1,N+1));
dist_trav_tot = Qs_nsc(jointi.pelvis.tx,end) -  Qs_nsc(jointi.pelvis.tx,1);
vel_aver_tot = dist_trav_tot/tf;
opti.subject_to(vel_aver_tot - S.v_tgt == 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale cost function
Jall_sc = sum(Jall)/dist_trav_tot;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create NLP solver
opti.minimize(Jall_sc);
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.mu_strategy           = 'adaptive';
options.ipopt.max_iter              = 10000;
options.ipopt.linear_solver         = S.linear_solver;
options.ipopt.tol                   = 1*10^(-S.tol_ipopt);
options.ipopt.constr_viol_tol       = 1*10^(-S.tol_ipopt);
opti.solver('ipopt', options);
% Create and save diary
OutFolder = fullfile(pathRepo,'Results',S.ResultsFolder);
if ~isfolder(OutFolder)
    mkdir(OutFolder);
end
Outname = fullfile(OutFolder,[S.savename '_log.txt']);
diary(Outname);
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
setup.tolerance.ipopt = S.tol_ipopt;
setup.bounds = bounds;
setup.scaling = scaling;
setup.guess = guess;

%% Save the results
% indicate that ExoVect is part of optimization results
if S.OptTexo_Ankle.Bool
    ExoVect = 'Optimization';
end

Outname = fullfile(OutFolder,[S.savename '.mat']);
Sopt = S;
save(Outname,'w_opt','stats','setup','Sopt','ExoVect');
end

