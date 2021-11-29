%%  Three-dimensional muscle-driven predictive simulations of human gaits

clear all; close all; clc;

% Name to save results
S.savename = 'TestPP_out';
S.loadname = 'TestPP';
S.ResultsFolder = 'Debug';

% flow control
S.Flow.solveProblem     = 0;   % set to 1 to solve problem
S.Flow.analyseResults   = 1;   % set to 1 to analyze results
S.Flow.loadResults      = 1;   % set to 1 to load results
S.Flow.saveResults      = 1;   % set to 1 to save sens. results
S.Flow.checkBoundsIG    = 0;   % set to 1 to visualize guess-bounds
S.Flow.writeIKmotion    = 1;   % set to 1 to write .mot file

% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.W.E       = 500;      % weight metabolic energy rate
S.W.Ak      = 50000;    % weight joint accelerations
S.W.ArmE    = 1000000;  % weight arm excitations
S.W.passMom = 1000;     % weight passive torques
S.W.A       = 2000;     % weight muscle activations
S.W.exp_E   = 2;        % power metabolic energy
S.W.Mtp     = 1000;     % weight mtp excitations
S.W.u       = 0.001;    % weight on excitations arms actuators
S.IGsel     = 2;        % initial guess identifier
S.ContactID = 1;        % contact model identifier
S.IGmodeID  = 1;        % initial guess mode identifier
S.IGcase    = 0;        % initial guess case identifier
S.h_weak    = 0;        % weakness hip actuators
S.Max_s     = 0;        % maximal contraction velocity identifier
S.pf_weak   = 0;        % weakness ankle plantaflexors
S.mE        = 0;        % metabolic energy model identifier
S.coCont    = 0;        % co-contraction identifier

% settings for exoskeleton control
S.ExoBool       = 1;
S.ExoScale      = 1;
S.DataSet       = 'PoggenSee2020_AFO';

% save name initial guess
S.savename_ig =[];


% ipopt options
S.linear_solver = 'mumps';
S.tol_ipopt     = 4;

% external function 
S.ExternalFunc = 'PredSim_mtp_cm1.dll';
S.ExternalFunc2 = 'PredSim_mtp_pp_cm1.dll';
%% User inputs (typical settings structure)

% flow control
solveProblem    = S.Flow.solveProblem; % set to 1 to solve problem
analyseResults  = S.Flow.analyseResults; % set to 1 to analyze results
loadResults     = S.Flow.loadResults; % set to 1 to load results
saveResults     = S.Flow.saveResults; % set to 1 to save sens. results
checkBoundsIG   = S.Flow.checkBoundsIG; % set to 1 to visualize guess-bounds
writeIKmotion   = S.Flow.writeIKmotion; % set to 1 to write .mot file

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
cm          = S.ContactID;  % contact model identifier
IGm         = S.IGmodeID;   % initial guess mode identifier
IGcase      = S.IGcase;     % initial guess case identifier
h_weak      = S.h_weak;     % weakness hip actuators
vMax_s      = S.Max_s;      % maximal contraction velocity identifier
pf_weak     = S.pf_weak;    % weakness ankle plantaflexors
mE          = S.mE;         % metabolic energy model identifier
coCont      = S.coCont;     % co-contraction identifier

% identifier for EMG load
savename_ig = S.savename_ig;

% ipopt options
tol_ipopt       = S.tol_ipopt;
linear_solver   = S.linear_solver;

%% Settings

import casadi.*
subject = 'subject1';
parallelMode = 'thread';
NThreads = 4;

%% Select settings
% Fixed parameter
W.u = 0.001;
% The filename used to save the results depends on the settings
v_tgt_id = round(v_tgt,2);
savename = S.savename;

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
F1 = external('F',S.ExternalFunc2);
cd(pathmain);


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
body_mass = 62;
body_weight = body_mass*9.81;

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
load([pathmusclemodel,'/MTparameters_',subject,'_mtp.mat']);
MTparameters_m = [MTparameters(:,musi),MTparameters(:,musi)];
% Indices of the muscles actuating the different joints for later use
pathpolynomial = [pathRepo,'/Polynomials'];
addpath(genpath(pathpolynomial));
tl = load([pathpolynomial,'/muscle_spanning_joint_INFO_',subject,'_mtp.mat']);
[~,mai] = MomentArmIndices(muscleNames(1:end-3),...
    tl.muscle_spanning_joint_INFO(1:end-3,:));

% By default, the tendon stiffness is 35 and the shift is 0.
aTendon = 35*ones(NMuscle,1);
shift = zeros(NMuscle,1);

% Adjust the maximal isometric force of the hip actuators if needed.
if h_weak ~= 0
    MTparameters_m(1,[mai(1).mus.l,mai(1).mus.r]) = ...
        MTparameters_m(1,[mai(1).mus.l,mai(1).mus.r])-...
        h_weak/100*MTparameters_m(1,[mai(1).mus.l,mai(1).mus.r]);
end

% Adjust the maximal isometric force of the ankle plantarflexors if needed.
if pf_weak ~= 0
    idx_FD = find(strcmp(muscleNames,'flex_dig_r'));
    idx_FH = find(strcmp(muscleNames,'flex_hal_r'));
    idx_GL = find(strcmp(muscleNames,'lat_gas_r'));
    idx_GM = find(strcmp(muscleNames,'med_gas_r'));
    idx_PB = find(strcmp(muscleNames,'per_brev_r'));
    idx_PL = find(strcmp(muscleNames,'per_long_r'));
    idx_SO = find(strcmp(muscleNames,'soleus_r'));
    idx_TP = find(strcmp(muscleNames,'tib_post_r'));
    idx_pf = [idx_FD,idx_FH,idx_GL,idx_GM,idx_PB,idx_PL,idx_SO,idx_TP];
    idx_pf_all = [idx_pf,idx_pf+NMuscle/2];
    MTparameters_m(1,idx_pf_all) = MTparameters_m(1,idx_pf_all)-...
        pf_weak/100*MTparameters_m(1,idx_pf_all);
end

% Adjust the maximum contraction velocities if needed.
if vMax_s == 1
    % Maximum contraction velocities * 2
    MTparameters_m(end,:) = MTparameters_m(end,:).*2;
elseif vMax_s == 2
    % Maximum contraction velocities * 1.5
    MTparameters_m(end,:) = MTparameters_m(end,:).*1.5;
end

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
addpath(genpath(pathCasADiFunctions));
% We load some variables for the polynomial approximations
load([pathpolynomial,'/muscle_spanning_joint_INFO_',subject,'_mtp.mat']);
load([pathpolynomial,'/MuscleInfo_',subject,'_mtp.mat']);
% For the polynomials, we want all independent muscles. So we do not need
% the muscles from both legs, since we assume bilateral symmetry, but want
% all muscles from the back (indices 47:49).
musi_pol = [musi,47,48,49];
NMuscle_pol = NMuscle/2+3;
CasADiFunctions_all_mtp

%% Passive joint torques
% We extract the parameters for the passive torques of the lower limbs and
% the trunk
pathPassiveMoments = [pathRepo,'/PassiveMoments'];
addpath(genpath(pathPassiveMoments));
PassiveMomentsData

stiffnessArm = 0;
dampingArm = 0.1;
stiffnessMtp = 1.5/(pi/180)/5;
dampingMtp = 0.5;

%% Experimental data
% We extract experimental data to set bounds and initial guesses if needed
pathData = [pathRepo,'/OpenSimModel/',subject];
joints = {'pelvis_tilt','pelvis_list','pelvis_rotation','pelvis_tx',...
    'pelvis_ty','pelvis_tz','hip_flexion_l','hip_adduction_l',...
    'hip_rotation_l','hip_flexion_r','hip_adduction_r','hip_rotation_r',...
    'knee_angle_l','knee_angle_r','ankle_angle_l','ankle_angle_r',...
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
    p = mfilename('fullpath');
    [~,namescript,~] = fileparts(p);
    pathIK = [pathRepo,'/Results/',namescript,'/IK',savename_ig,'.mot'];
    Qs_ig = getIK(pathIK,joints);
    % When saving the results, we save a full gait cycle (2*N) so here we
    % only select 1:N to have half a gait cycle
    Qs_ig_sel.allfilt = Qs_ig.allfilt(1:N,:);
    Qs_ig_sel.time = Qs_ig.time(1:N,:);
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
    guess = getGuess_QR_opti_int(N,nq,NMuscle,scaling,v_tgt,jointi,d);
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
        guess = ...
            getGuess_DI_t(Qs_ig_sel,nq,N,time_IC,NMuscle,jointi,scaling);
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
% This allows visualizing the initial guess and the bounds
if checkBoundsIG
    pathPlots = [pathRepo,'/Plots'];
    addpath(genpath(pathPlots));
    plot_BoundsVSInitialGuess_all_mtp
end

%% exoskeleton torques
ExoControl = [];
body_mass = 62;
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


%% Formulate the NLP
if solveProblem
    % Create opti instance
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
    % Parallel formulation
    % Define CasADi variables for static parameters
    tfk = MX.sym('tfk');
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
            lMTj_lr,vMTj_lr,tensions,aTendon,shift);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get metabolic energy rate if in the cost function        
        % Get muscle fiber lengths
        [~,lMtildej] = f_FiberLength_TendonForce_tendon(...
            FTtildekj_nsc(:,j+1),lMTj_lr,aTendon,shift);
        % Get muscle fiber velocities
        [vMj,~] = f_FiberVelocity_TendonForce_tendon(...
            FTtildekj_nsc(:,j+1),dFTtildej_nsc(:,j),...
            lMTj_lr,vMTj_lr,aTendon,shift);
        % Get metabolic energy rate Bhargava et al. (2004)
            [e_totj,~,~,~,~,~] = fgetMetabolicEnergySmooth2004all(...
                akj(:,j+1),akj(:,j+1),lMtildej,...
                vMj,Fcej,Fpassj,massMj,pctsts,Fisoj,...
                MTparameters_m(1,:)',body_mass,10);
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get passive joint torques
        Tau_passj.hip.flex.l    = f_PassiveMoments(k_pass.hip.flex,...
            theta.pass.hip.flex,Qskj_nsc(jointi.hip_flex.l,j+1),...
            Qdotskj_nsc(jointi.hip_flex.l,j+1));
        Tau_passj.hip.flex.r    = f_PassiveMoments(k_pass.hip.flex,...
            theta.pass.hip.flex,Qskj_nsc(jointi.hip_flex.r,j+1),...
            Qdotskj_nsc(jointi.hip_flex.r,j+1));
        Tau_passj.hip.add.l     = f_PassiveMoments(k_pass.hip.add,...
            theta.pass.hip.add,Qskj_nsc(jointi.hip_add.l,j+1),...
            Qdotskj_nsc(jointi.hip_add.l,j+1));
        Tau_passj.hip.add.r     = f_PassiveMoments(k_pass.hip.add,...
            theta.pass.hip.add,Qskj_nsc(jointi.hip_add.r,j+1),...
            Qdotskj_nsc(jointi.hip_add.r,j+1));
        Tau_passj.hip.rot.l     = f_PassiveMoments(k_pass.hip.rot,...
            theta.pass.hip.rot,Qskj_nsc(jointi.hip_rot.l,j+1),...
            Qdotskj_nsc(jointi.hip_rot.l,j+1));
        Tau_passj.hip.rot.r     = f_PassiveMoments(k_pass.hip.rot,...
            theta.pass.hip.rot,Qskj_nsc(jointi.hip_rot.r,j+1),...
            Qdotskj_nsc(jointi.hip_rot.r,j+1));
        Tau_passj.knee.l        = f_PassiveMoments(k_pass.knee,...
            theta.pass.knee,Qskj_nsc(jointi.knee.l,j+1),...
            Qdotskj_nsc(jointi.knee.l,j+1));
        Tau_passj.knee.r        = f_PassiveMoments(k_pass.knee,...
            theta.pass.knee,Qskj_nsc(jointi.knee.r,j+1),...
            Qdotskj_nsc(jointi.knee.r,j+1));
        Tau_passj.ankle.l       = f_PassiveMoments(k_pass.ankle,...
            theta.pass.ankle,Qskj_nsc(jointi.ankle.l,j+1),...
            Qdotskj_nsc(jointi.ankle.l,j+1));
        Tau_passj.ankle.r       = f_PassiveMoments(k_pass.ankle,...
            theta.pass.ankle,Qskj_nsc(jointi.ankle.r,j+1),...
            Qdotskj_nsc(jointi.ankle.r,j+1));
        Tau_passj.subt.l       = f_PassiveMoments(k_pass.subt,...
            theta.pass.subt,Qskj_nsc(jointi.subt.l,j+1),...
            Qdotskj_nsc(jointi.subt.l,j+1));
        Tau_passj.subt.r       = f_PassiveMoments(k_pass.subt,...
            theta.pass.subt,Qskj_nsc(jointi.subt.r,j+1),...
            Qdotskj_nsc(jointi.subt.r,j+1));
        Tau_passj.trunk.ext     = f_PassiveMoments(k_pass.trunk.ext,...
            theta.pass.trunk.ext,Qskj_nsc(jointi.trunk.ext,j+1),...
            Qdotskj_nsc(jointi.trunk.ext,j+1));
        Tau_passj.trunk.ben     = f_PassiveMoments(k_pass.trunk.ben,...
            theta.pass.trunk.ben,Qskj_nsc(jointi.trunk.ben,j+1),...
            Qdotskj_nsc(jointi.trunk.ben,j+1));
        Tau_passj.trunk.rot     = f_PassiveMoments(k_pass.trunk.rot,...
            theta.pass.trunk.rot,Qskj_nsc(jointi.trunk.rot,j+1),...
            Qdotskj_nsc(jointi.trunk.rot,j+1));
        
        Tau_passj.sh_flex.l = f_passiveTATorques(stiffnessArm, dampingArm, ...
            Qskj_nsc(jointi.sh_flex.l,j+1), Qdotskj_nsc(jointi.sh_flex.l,j+1));
        Tau_passj.sh_add.l = f_passiveTATorques(stiffnessArm, dampingArm, ...
            Qskj_nsc(jointi.sh_add.l,j+1), Qdotskj_nsc(jointi.sh_add.l,j+1));
        Tau_passj.sh_rot.l = f_passiveTATorques(stiffnessArm, dampingArm, ...
            Qskj_nsc(jointi.sh_rot.l,j+1), Qdotskj_nsc(jointi.sh_rot.l,j+1));
        Tau_passj.sh_flex.r = f_passiveTATorques(stiffnessArm, dampingArm, ...
            Qskj_nsc(jointi.sh_flex.r,j+1), Qdotskj_nsc(jointi.sh_flex.r,j+1));
        Tau_passj.sh_add.r = f_passiveTATorques(stiffnessArm, dampingArm, ...
            Qskj_nsc(jointi.sh_add.r,j+1), Qdotskj_nsc(jointi.sh_add.r,j+1));
        Tau_passj.sh_rot.r = f_passiveTATorques(stiffnessArm, dampingArm, ...
            Qskj_nsc(jointi.sh_rot.r,j+1), Qdotskj_nsc(jointi.sh_rot.r,j+1));
        Tau_passj.elb.l = f_passiveTATorques(stiffnessArm, dampingArm, ...
            Qskj_nsc(jointi.elb.l,j+1), Qdotskj_nsc(jointi.elb.l,j+1));
        Tau_passj.elb.r = f_passiveTATorques(stiffnessArm, dampingArm, ...
            Qskj_nsc(jointi.elb.r,j+1), Qdotskj_nsc(jointi.elb.r,j+1));
        Tau_passj.arm = [Tau_passj.sh_flex.l, Tau_passj.sh_add.l, ...
            Tau_passj.sh_rot.l, Tau_passj.sh_flex.r, Tau_passj.sh_add.r, ...
            Tau_passj.sh_rot.r, Tau_passj.elb.l, Tau_passj.elb.r];
        
        Tau_passj.mtp.l = f_passiveTATorques(stiffnessMtp, dampingMtp, ...
            Qskj_nsc(jointi.mtp.l,j+1), Qdotskj_nsc(jointi.mtp.l,j+1));
        Tau_passj.mtp.r = f_passiveTATorques(stiffnessMtp, dampingMtp, ...
            Qskj_nsc(jointi.mtp.r,j+1), Qdotskj_nsc(jointi.mtp.r,j+1));
        Tau_passj.mtp.all = [Tau_passj.mtp.l, Tau_passj.mtp.r];
        
        Tau_passj_all = [Tau_passj.hip.flex.l,Tau_passj.hip.flex.r,...
            Tau_passj.hip.add.l,Tau_passj.hip.add.r,...
            Tau_passj.hip.rot.l,Tau_passj.hip.rot.r,...
            Tau_passj.knee.l,Tau_passj.knee.r,Tau_passj.ankle.l,...
            Tau_passj.ankle.r,Tau_passj.subt.l,Tau_passj.subt.r,...
            Tau_passj.mtp.all,Tau_passj.trunk.ext,Tau_passj.trunk.ben,...
            Tau_passj.trunk.rot,Tau_passj.arm]';
        
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
        if W.E == 0
            J = J + 1*(...
                W.A*B(j+1)      *(f_J92(akj(:,j+1)'))*h + ...
                W.ArmE*B(j+1)   *(f_J8(e_ak))*h +...
                W.Ak*B(j+1)     *(f_J21(Aj(residuals_noarmsi,j)))*h + ...
                W.passMom*B(j+1)*(f_J15(Tau_passj_all))*h + ...
                W.u*B(j+1)      *(f_J92(vAk))*h + ...
                W.u*B(j+1)      *(f_J92(dFTtildej(:,j)))*h + ...
                W.u*B(j+1)      *(f_J8(Aj(armsi,j)))*h);
        elseif W.A == 0
            J = J + 1*(...
                W.E*B(j+1)      *(f_J92exp(e_totj,exp_E))/body_mass*h + ...
                W.ArmE*B(j+1)   *(f_J8(e_ak))*h +...
                W.Ak*B(j+1)     *(f_J21(Aj(residuals_noarmsi,j)))*h + ...
                W.passMom*B(j+1)*(f_J15(Tau_passj_all))*h + ...
                W.u*B(j+1)      *(f_J92(vAk))*h + ...
                W.u*B(j+1)      *(f_J92(dFTtildej(:,j)))*h + ...
                W.u*B(j+1)      *(f_J8(Aj(armsi,j)))*h);
        elseif W.passMom == 0
            J = J + 1*(...
                W.E*B(j+1)      *(f_J92exp(e_totj,exp_E))/body_mass*h + ...
                W.A*B(j+1)      *(f_J92(akj(:,j+1)'))*h + ...
                W.ArmE*B(j+1)   *(f_J8(e_ak))*h +...
                W.Ak*B(j+1)     *(f_J21(Aj(residuals_noarmsi,j)))*h + ...
                W.u*B(j+1)      *(f_J92(vAk))*h + ...
                W.u*B(j+1)      *(f_J92(dFTtildej(:,j)))*h + ...
                W.u*B(j+1)      *(f_J8(Aj(armsi,j)))*h);
        else
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
        end
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
        eq_constr{end+1} = Tj(jointi.hip_flex.l,1)-(T_hip_flex_l + ...
            Tau_passj.hip.flex.l);
        % Hip flexion, right
        Ft_hip_flex_r   = FTj(mai(1).mus.r',1);
        T_hip_flex_r    = f_T27(MAj.hip_flex.r,Ft_hip_flex_r);
        eq_constr{end+1} = Tj(jointi.hip_flex.r,1)-(T_hip_flex_r + ...
            Tau_passj.hip.flex.r);
        % Hip adduction, left
        Ft_hip_add_l    = FTj(mai(2).mus.l',1);
        T_hip_add_l     = f_T27(MAj.hip_add.l,Ft_hip_add_l);
        eq_constr{end+1} = Tj(jointi.hip_add.l,1)-(T_hip_add_l + ...
            Tau_passj.hip.add.l);
        % Hip adduction, right
        Ft_hip_add_r    = FTj(mai(2).mus.r',1);
        T_hip_add_r     = f_T27(MAj.hip_add.r,Ft_hip_add_r);
        eq_constr{end+1} = Tj(jointi.hip_add.r,1)-(T_hip_add_r + ...
            Tau_passj.hip.add.r);
        % Hip rotation, left
        Ft_hip_rot_l    = FTj(mai(3).mus.l',1);
        T_hip_rot_l     = f_T27(MAj.hip_rot.l,Ft_hip_rot_l);
        eq_constr{end+1} = Tj(jointi.hip_rot.l,1)-(T_hip_rot_l + ...
            Tau_passj.hip.rot.l);
        % Hip rotation, right
        Ft_hip_rot_r    = FTj(mai(3).mus.r',1);
        T_hip_rot_r     = f_T27(MAj.hip_rot.r,Ft_hip_rot_r);
        eq_constr{end+1} = Tj(jointi.hip_rot.r,1)-(T_hip_rot_r + ...
            Tau_passj.hip.rot.r);
        % Knee, left
        Ft_knee_l       = FTj(mai(4).mus.l',1);
        T_knee_l        = f_T13(MAj.knee.l,Ft_knee_l);
        eq_constr{end+1} = Tj(jointi.knee.l,1)-(T_knee_l + ...
            Tau_passj.knee.l);
        % Knee, right
        Ft_knee_r       = FTj(mai(4).mus.r',1);
        T_knee_r        = f_T13(MAj.knee.r,Ft_knee_r);
        eq_constr{end+1} = Tj(jointi.knee.r,1)-(T_knee_r + ...
            Tau_passj.knee.r);
        % Ankle, left
        Ft_ankle_l      = FTj(mai(5).mus.l',1);
        T_ankle_l       = f_T12(MAj.ankle.l,Ft_ankle_l);
        eq_constr{end+1} = Tj(jointi.ankle.l,1)-(T_ankle_l + ...
            Tau_passj.ankle.l );
        % Ankle, right
        Ft_ankle_r      = FTj(mai(5).mus.r',1);
        T_ankle_r       = f_T12(MAj.ankle.r,Ft_ankle_r);
        eq_constr{end+1} = Tj(jointi.ankle.r,1)-(T_ankle_r + ...
            Tau_passj.ankle.r);
        % Subtalar, left
        Ft_subt_l       = FTj(mai(6).mus.l',1);
        T_subt_l        = f_T12(MAj.subt.l,Ft_subt_l);
        eq_constr{end+1} = Tj(jointi.subt.l,1)-(T_subt_l + ...
            Tau_passj.subt.l + Texok(1));
        % Subtalar, right
        Ft_subt_r       = FTj(mai(6).mus.r',1);
        T_subt_r        = f_T12(MAj.subt.r,Ft_subt_r);
        eq_constr{end+1} = Tj(jointi.subt.r,1)-(T_subt_r + ...
            Tau_passj.subt.r + Texok(2));
        % Lumbar extension
        Ft_trunk_ext    = FTj([mai(8).mus.l,mai(8).mus.r]',1);
        T_trunk_ext     = f_T6(MAj.trunk_ext,Ft_trunk_ext);
        eq_constr{end+1} = Tj(jointi.trunk.ext,1)-(T_trunk_ext + ...
            Tau_passj.trunk.ext);
        % Lumbar bending
        Ft_trunk_ben    = FTj([mai(9).mus.l,mai(9).mus.r]',1);
        T_trunk_ben     = f_T6(MAj.trunk_ben,Ft_trunk_ben);
        eq_constr{end+1} = Tj(jointi.trunk.ben,1)-(T_trunk_ben + ...
            Tau_passj.trunk.ben);
        % Lumbar rotation
        Ft_trunk_rot    = FTj([mai(10).mus.l,mai(10).mus.r]',1);
        T_trunk_rot     = f_T6(MAj.trunk_rot,Ft_trunk_rot);
        eq_constr{end+1} = Tj(jointi.trunk.rot,1)-(T_trunk_rot + ...
            Tau_passj.trunk.rot);
        % Torque-driven joint torques for the arms
        % Arms
        eq_constr{end+1} = Tj(armsi,1)/scaling.ArmTau - (a_akj(:,j+1) + ...
            (Tau_passj.arm)'/scaling.ArmTau);
        % Mtp
        eq_constr{end+1} = Tj(mtpi,1)/scaling.MtpTau - (a_mtpkj(:,j+1) + ...
            (Tau_passj.mtp.all)'/scaling.MtpTau);
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
    p = mfilename('fullpath');
    [~,namescript,~] = fileparts(p);
    pathresults = [pathRepo,'/Results'];
    if ~(exist([pathresults,'/',namescript],'dir')==7)
        mkdir(pathresults,namescript);
    end
    if (exist([pathresults,'/',namescript,'/D',savename],'file')==2)
        delete ([pathresults,'/',namescript,'/D',savename])
    end
    diary([pathresults,'/',namescript,'/D',savename]);
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

    % save the results in the right folder:
    OutFolder = fullfile(pathRepo,'Results',S.ResultsFolder);
    if ~isfolder(OutFolder)
        mkdir(OutFolder);
    end
    Outname = fullfile(OutFolder,[S.savename '.mat']); 
    Sopt = S;
    save(Outname,'w_opt','stats','setup','Sopt','ExoControl');  

    % Save results and set
    %save([pathresults,'/',namescript,'/w',savename],'w_opt');
    %save([pathresults,'/',namescript,'/stats',savename],'stats');
end

%% Analyze results
if analyseResults
    %% Load results
    if loadResults
        OutFolder = fullfile(pathRepo,'Results',S.ResultsFolder);
        Outname = fullfile(OutFolder,[S.loadname '.mat']); 
        load(Outname);
    end
    NParameters = 1;
    tf_opt = w_opt(1:NParameters);
    starti = NParameters+1;
    a_opt = reshape(w_opt(starti:starti+NMuscle*(N+1)-1),NMuscle,N+1)';
    starti = starti + NMuscle*(N+1);
    a_col_opt = reshape(w_opt(starti:starti+NMuscle*(d*N)-1),NMuscle,d*N)';
    starti = starti + NMuscle*(d*N);
    FTtilde_opt = reshape(w_opt(starti:starti+NMuscle*(N+1)-1),NMuscle,N+1)';
    starti = starti + NMuscle*(N+1);
    FTtilde_col_opt =reshape(w_opt(starti:starti+NMuscle*(d*N)-1),NMuscle,d*N)';
    starti = starti + NMuscle*(d*N);
    Qs_opt = reshape(w_opt(starti:starti+nq.all*(N+1)-1),nq.all,N+1)';
    starti = starti + nq.all*(N+1);
    Qs_col_opt = reshape(w_opt(starti:starti+nq.all*(d*N)-1),nq.all,d*N)';
    starti = starti + nq.all*(d*N);
    Qdots_opt = reshape(w_opt(starti:starti+nq.all*(N+1)-1),nq.all,N+1)';
    starti = starti + nq.all*(N+1);
    Qdots_col_opt = reshape(w_opt(starti:starti+nq.all*(d*N)-1),nq.all,d*N)';
    starti = starti + nq.all*(d*N);
    a_a_opt = reshape(w_opt(starti:starti+nq.arms*(N+1)-1),nq.arms,N+1)';
    starti = starti + nq.arms*(N+1);
    a_a_col_opt = reshape(w_opt(starti:starti+nq.arms*(d*N)-1),nq.arms,d*N)';
    starti = starti + nq.arms*(d*N);
    a_mtp_opt = reshape(w_opt(starti:starti+nq.mtp*(N+1)-1),nq.mtp,N+1)';
    starti = starti + nq.mtp*(N+1);
    a_mtp_col_opt = reshape(w_opt(starti:starti+nq.mtp*(d*N)-1),nq.mtp,d*N)';
    starti = starti + nq.mtp*(d*N);
    vA_opt = reshape(w_opt(starti:starti+NMuscle*N-1),NMuscle,N)';
    starti = starti + NMuscle*N;
    e_a_opt = reshape(w_opt(starti:starti+nq.arms*N-1),nq.arms,N)';
    starti = starti + nq.arms*N;
    e_mtp_opt = reshape(w_opt(starti:starti+nq.mtp*N-1),nq.mtp,N)';
    starti = starti + nq.mtp*N;
    dFTtilde_col_opt=reshape(w_opt(starti:starti+NMuscle*(d*N)-1),NMuscle,d*N)';
    starti = starti + NMuscle*(d*N);
    qdotdot_col_opt =reshape(w_opt(starti:starti+nq.all*(d*N)-1),nq.all,(d*N))';
    starti = starti + nq.all*(d*N);
    if starti - 1 ~= length(w_opt)
        disp('error when extracting results')
    end
    % Combine results at mesh and collocation points
    a_mesh_col_opt=zeros(N*(d+1)+1,NMuscle);
    a_mesh_col_opt(1:(d+1):end,:)= a_opt;
    FTtilde_mesh_col_opt=zeros(N*(d+1)+1,NMuscle);
    FTtilde_mesh_col_opt(1:(d+1):end,:)= FTtilde_opt;
    Qs_mesh_col_opt=zeros(N*(d+1)+1,nq.all);
    Qs_mesh_col_opt(1:(d+1):end,:)= Qs_opt;
    Qdots_mesh_col_opt=zeros(N*(d+1)+1,nq.all);
    Qdots_mesh_col_opt(1:(d+1):end,:)= Qdots_opt;
    a_a_mesh_col_opt=zeros(N*(d+1)+1,nq.arms);
    a_a_mesh_col_opt(1:(d+1):end,:)= a_a_opt;
    a_mtp_mesh_col_opt=zeros(N*(d+1)+1,nq.mtp);
    a_mtp_mesh_col_opt(1:(d+1):end,:)= a_mtp_opt;
    for k=1:N
        rangei = k*(d+1)-(d-1):k*(d+1);
        rangebi = (k-1)*d+1:k*d;
        a_mesh_col_opt(rangei,:) = a_col_opt(rangebi,:);
        FTtilde_mesh_col_opt(rangei,:) = FTtilde_col_opt(rangebi,:);
        Qs_mesh_col_opt(rangei,:) = Qs_col_opt(rangebi,:);
        Qdots_mesh_col_opt(rangei,:) = Qdots_col_opt(rangebi,:);
        a_a_mesh_col_opt(rangei,:) = a_a_col_opt(rangebi,:);
        a_mtp_mesh_col_opt(rangei,:) = a_mtp_col_opt(rangebi,:);
    end
    
    %% Unscale results
    % States at mesh points
    % Qs (1:N-1)
    q_opt_unsc.rad = Qs_opt(1:end-1,:).*repmat(...
        scaling.Qs,size(Qs_opt(1:end-1,:),1),1);
    % Convert in degrees
    q_opt_unsc.deg = q_opt_unsc.rad;
    q_opt_unsc.deg(:,roti) = q_opt_unsc.deg(:,roti).*180/pi;
    % Qs (1:N)
    q_opt_unsc_all.rad = Qs_opt.*repmat(scaling.Qs,size(Qs_opt,1),1);
    % Convert in degrees
    q_opt_unsc_all.deg = q_opt_unsc_all.rad;
    q_opt_unsc_all.deg(:,roti) = q_opt_unsc_all.deg(:,roti).*180/pi;
    % Qdots (1:N-1)
    qdot_opt_unsc.rad = Qdots_opt(1:end-1,:).*repmat(...
        scaling.Qdots,size(Qdots_opt(1:end-1,:),1),1);
    % Convert in degrees
    qdot_opt_unsc.deg = qdot_opt_unsc.rad;
    qdot_opt_unsc.deg(:,roti) = qdot_opt_unsc.deg(:,roti).*180/pi;
    % Qdots (1:N)
    qdot_opt_unsc_all.rad =Qdots_opt.*repmat(scaling.Qdots,size(Qdots_opt,1),1);
    % Muscle activations (1:N-1)
    a_opt_unsc = a_opt(1:end-1,:).*repmat(...
        scaling.a,size(a_opt(1:end-1,:),1),size(a_opt,2));
    % Muscle-tendon forces (1:N-1)
    FTtilde_opt_unsc = FTtilde_opt(1:end-1,:).*repmat(...
        scaling.FTtilde,size(FTtilde_opt(1:end-1,:),1),1);
    % Arm activations (1:N-1)
    a_a_opt_unsc = a_a_opt(1:end-1,:);
    % Arm activations (1:N)
    a_a_opt_unsc_all = a_a_opt;
    % Mtp activations (1:N-1)
    a_mtp_opt_unsc = a_mtp_opt(1:end-1,:);
    % Mtp activations (1:N)
    a_mtp_opt_unsc_all = a_mtp_opt;
    % Controls at mesh points
    % Time derivative of muscle activations (states)
    vA_opt_unsc = vA_opt.*repmat(scaling.vA,size(vA_opt,1),size(vA_opt,2));
    % Get muscle excitations from time derivative of muscle activations
    e_opt_unsc = computeExcitationRaasch(a_opt_unsc,vA_opt_unsc,...
        ones(1,NMuscle)*tdeact,ones(1,NMuscle)*tact);
    % Arm excitations
    e_a_opt_unsc = e_a_opt;
    % Mtp excitations
    e_mtp_opt_unsc = e_mtp_opt;
    % States at collocation points
    % Qs
    q_col_opt_unsc.rad = Qs_col_opt.*repmat(scaling.Qs,size(Qs_col_opt,1),1);
    % Convert in degrees
    q_col_opt_unsc.deg = q_col_opt_unsc.rad;
    q_col_opt_unsc.deg(:,roti) = q_col_opt_unsc.deg(:,roti).*180/pi;
    % Qdots
    qdot_col_opt_unsc.rad = Qdots_col_opt.*repmat(...
        scaling.Qdots,size(Qdots_col_opt,1),1);
    % Convert in degrees
    qdot_col_opt_unsc.deg = qdot_col_opt_unsc.rad;
    qdot_col_opt_unsc.deg(:,roti) = qdot_col_opt_unsc.deg(:,roti).*180/pi;
    % Muscle activations
    a_col_opt_unsc = a_col_opt.*repmat(...
        scaling.a,size(a_col_opt,1),size(a_col_opt,2));
    % Muscle-tendon forces
    FTtilde_col_opt_unsc = FTtilde_col_opt.*repmat(...
        scaling.FTtilde,size(FTtilde_col_opt,1),1);
    % Arm activations
    a_a_col_opt_unsc = a_a_col_opt;
    % Mtp activations
    a_mtp_col_opt_unsc = a_mtp_col_opt;
    % "Slack" controls at collocation points
    % Time derivative of Qdots
    qdotdot_col_opt_unsc.rad = ...
        qdotdot_col_opt.*repmat(scaling.Qdotdots,size(qdotdot_col_opt,1),1);
    % Convert in degrees
    qdotdot_col_opt_unsc.deg = qdotdot_col_opt_unsc.rad;
    qdotdot_col_opt_unsc.deg(:,roti) = qdotdot_col_opt_unsc.deg(:,roti).*180/pi;
    % Time derivative of muscle-tendon forces
    dFTtilde_col_opt_unsc = dFTtilde_col_opt.*repmat(...
        scaling.dFTtilde,size(dFTtilde_col_opt,1),size(dFTtilde_col_opt,2));
    dFTtilde_opt_unsc = dFTtilde_col_opt_unsc(d:d:end,:);
    
    %% Time grid
    % Mesh points
    tgrid = linspace(0,tf_opt,N+1);
    dtime = zeros(1,d+1);
    for i=1:4
        dtime(i)=tau_root(i)*(tf_opt/N);
    end
    % Mesh points and collocation points
    tgrid_ext = zeros(1,(d+1)*N+1);
    for i=1:N
        tgrid_ext(((i-1)*4+1):1:i*4)=tgrid(i)+dtime;
    end
    tgrid_ext(end)=tf_opt;
    
    %% Joint torques and ground reaction forces at mesh points (N-1), except #1
    Xk_Qs_Qdots_opt             = zeros(N,2*nq.all);
    Xk_Qs_Qdots_opt(:,1:2:end)  = q_opt_unsc_all.rad(2:end,:);
    Xk_Qs_Qdots_opt(:,2:2:end)  = qdot_opt_unsc_all.rad(2:end,:);
    Xk_Qdotdots_opt             = qdotdot_col_opt_unsc.rad(d:d:end,:);
    Foutk_opt = zeros(N,nq.all+NGRF+NcalcOrall);
    Tau_passk_opt_all = zeros(N,nq.all-nq.abs);
    for i = 1:N
        [res] = F1([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)']);
        Foutk_opt(i,:) = full(res);
        
        Tau_passk_opt.hip.flex.l = full(f_PassiveMoments(k_pass.hip.flex,...
            theta.pass.hip.flex,q_opt_unsc_all.rad(i+1,jointi.hip_flex.l),...
            qdot_opt_unsc_all.rad(i+1,jointi.hip_flex.l)));
        Tau_passk_opt.hip.flex.r = full(f_PassiveMoments(k_pass.hip.flex,...
            theta.pass.hip.flex,q_opt_unsc_all.rad(i+1,jointi.hip_flex.r),...
            qdot_opt_unsc_all.rad(i+1,jointi.hip_flex.r)));
        Tau_passk_opt.hip.add.l = full(f_PassiveMoments(k_pass.hip.add,...
            theta.pass.hip.add,q_opt_unsc_all.rad(i+1,jointi.hip_add.l),...
            qdot_opt_unsc_all.rad(i+1,jointi.hip_add.l)));
        Tau_passk_opt.hip.add.r = full(f_PassiveMoments(k_pass.hip.add,...
            theta.pass.hip.add,q_opt_unsc_all.rad(i+1,jointi.hip_add.r),...
            qdot_opt_unsc_all.rad(i+1,jointi.hip_add.r)));
        Tau_passk_opt.hip.rot.l = full(f_PassiveMoments(k_pass.hip.rot,...
            theta.pass.hip.rot,q_opt_unsc_all.rad(i+1,jointi.hip_rot.l),...
            qdot_opt_unsc_all.rad(i+1,jointi.hip_rot.l)));
        Tau_passk_opt.hip.rot.r = full(f_PassiveMoments(k_pass.hip.rot,...
            theta.pass.hip.rot,q_opt_unsc_all.rad(i+1,jointi.hip_rot.r),...
            qdot_opt_unsc_all.rad(i+1,jointi.hip_rot.r)));
        Tau_passk_opt.knee.l = full(f_PassiveMoments(k_pass.knee,...
            theta.pass.knee,q_opt_unsc_all.rad(i+1,jointi.knee.l),...
            qdot_opt_unsc_all.rad(i+1,jointi.knee.l)));
        Tau_passk_opt.knee.r = full(f_PassiveMoments(k_pass.knee,...
            theta.pass.knee,q_opt_unsc_all.rad(i+1,jointi.knee.r),...
            qdot_opt_unsc_all.rad(i+1,jointi.knee.r)));
        Tau_passk_opt.ankle.l = full(f_PassiveMoments(k_pass.ankle,...
            theta.pass.ankle,q_opt_unsc_all.rad(i+1,jointi.ankle.l),...
            qdot_opt_unsc_all.rad(i+1,jointi.ankle.l)));
        Tau_passk_opt.ankle.r = full(f_PassiveMoments(k_pass.ankle,...
            theta.pass.ankle,q_opt_unsc_all.rad(i+1,jointi.ankle.r),...
            qdot_opt_unsc_all.rad(i+1,jointi.ankle.r)));
        Tau_passk_opt.subt.l = full(f_PassiveMoments(k_pass.subt,...
            theta.pass.subt,q_opt_unsc_all.rad(i+1,jointi.subt.l),...
            qdot_opt_unsc_all.rad(i+1,jointi.subt.l)));
        Tau_passk_opt.subt.r = full(f_PassiveMoments(k_pass.subt,...
            theta.pass.subt,q_opt_unsc_all.rad(i+1,jointi.subt.r),...
            qdot_opt_unsc_all.rad(i+1,jointi.subt.r)));
        Tau_passk_opt.trunk.ext = full(f_PassiveMoments(k_pass.trunk.ext,...
            theta.pass.trunk.ext,q_opt_unsc_all.rad(i+1,jointi.trunk.ext),...
            qdot_opt_unsc_all.rad(i+1,jointi.trunk.ext)));
        Tau_passk_opt.trunk.ben = full(f_PassiveMoments(k_pass.trunk.ben,...
            theta.pass.trunk.ben,q_opt_unsc_all.rad(i+1,jointi.trunk.ben),...
            qdot_opt_unsc_all.rad(i+1,jointi.trunk.ben)));
        Tau_passk_opt.trunk.rot = full(f_PassiveMoments(k_pass.trunk.rot,...
            theta.pass.trunk.rot,q_opt_unsc_all.rad(i+1,jointi.trunk.rot),...
            qdot_opt_unsc_all.rad(i+1,jointi.trunk.rot)));
        
        Tau_passk_opt.sh_flex.l = full(f_passiveTATorques(stiffnessArm, dampingArm, ...
            q_opt_unsc_all.rad(i+1,jointi.sh_flex.l), qdot_opt_unsc_all.rad(i+1,jointi.sh_flex.l)));
        Tau_passk_opt.sh_add.l = full(f_passiveTATorques(stiffnessArm, dampingArm, ...
            q_opt_unsc_all.rad(i+1,jointi.sh_add.l), qdot_opt_unsc_all.rad(i+1,jointi.sh_add.l)));
        Tau_passk_opt.sh_rot.l = full(f_passiveTATorques(stiffnessArm, dampingArm, ...
            q_opt_unsc_all.rad(i+1,jointi.sh_rot.l), qdot_opt_unsc_all.rad(i+1,jointi.sh_rot.l)));
        Tau_passk_opt.sh_flex.r = full(f_passiveTATorques(stiffnessArm, dampingArm, ...
            q_opt_unsc_all.rad(i+1,jointi.sh_flex.r), qdot_opt_unsc_all.rad(i+1,jointi.sh_flex.r)));
        Tau_passk_opt.sh_add.r = full(f_passiveTATorques(stiffnessArm, dampingArm, ...
            q_opt_unsc_all.rad(i+1,jointi.sh_add.r), qdot_opt_unsc_all.rad(i+1,jointi.sh_add.r)));
        Tau_passk_opt.sh_rot.r = full(f_passiveTATorques(stiffnessArm, dampingArm, ...
            q_opt_unsc_all.rad(i+1,jointi.sh_rot.r), qdot_opt_unsc_all.rad(i+1,jointi.sh_rot.r)));
        Tau_passk_opt.elb.l = full(f_passiveTATorques(stiffnessArm, dampingArm, ...
            q_opt_unsc_all.rad(i+1,jointi.elb.l), qdot_opt_unsc_all.rad(i+1,jointi.elb.l)));
        Tau_passk_opt.elb.r = full(f_passiveTATorques(stiffnessArm, dampingArm, ...
            q_opt_unsc_all.rad(i+1,jointi.elb.r), qdot_opt_unsc_all.rad(i+1,jointi.elb.r)));
        Tau_passk_opt.arm(i,:) = [Tau_passk_opt.sh_flex.l, Tau_passk_opt.sh_add.l, ...
            Tau_passk_opt.sh_rot.l, Tau_passk_opt.sh_flex.r, Tau_passk_opt.sh_add.r, ...
            Tau_passk_opt.sh_rot.r, Tau_passk_opt.elb.l, Tau_passk_opt.elb.r];
        
        Tau_passk_opt.mtp.l = full(f_passiveTATorques(stiffnessMtp, dampingMtp, ...
            q_opt_unsc_all.rad(i+1,jointi.mtp.l), qdot_opt_unsc_all.rad(i+1,jointi.mtp.l)));
        Tau_passk_opt.mtp.r = full(f_passiveTATorques(stiffnessMtp, dampingMtp, ...
            q_opt_unsc_all.rad(i+1,jointi.mtp.r), qdot_opt_unsc_all.rad(i+1,jointi.mtp.r)));
        Tau_passk_opt.mtp.all(i,:) = [Tau_passk_opt.mtp.l, Tau_passk_opt.mtp.r];
        
        Tau_passk_opt_all(i,:) = [Tau_passk_opt.hip.flex.l,...
            Tau_passk_opt.hip.add.l,Tau_passk_opt.hip.rot.l,...
            Tau_passk_opt.hip.flex.r,Tau_passk_opt.hip.add.r,...
            Tau_passk_opt.hip.rot.r,Tau_passk_opt.knee.l,...
            Tau_passk_opt.knee.r,Tau_passk_opt.ankle.l,...
            Tau_passk_opt.ankle.r,Tau_passk_opt.subt.l,...
            Tau_passk_opt.subt.r,Tau_passk_opt.mtp.all(i,:),...
            Tau_passk_opt.trunk.ext,Tau_passk_opt.trunk.ben,...
            Tau_passk_opt.trunk.rot,Tau_passk_opt.arm(i,:),];
        
    end
    GRFk_opt = Foutk_opt(:,GRFi.all);
    assertArmTmaxk = max(max(abs(Foutk_opt(:,armsi)/scaling.ArmTau - ...
        (a_a_opt_unsc_all(2:end,:) + Tau_passk_opt.arm/scaling.ArmTau))));
    assertMtpTmaxk = max(max(abs(Foutk_opt(:,mtpi)/scaling.MtpTau - ...
        (a_mtp_opt_unsc_all(2:end,:) + Tau_passk_opt.mtp.all/scaling.MtpTau))));
    if assertArmTmaxk > 1*10^(-tol_ipopt)
        disp(['Issue when reconstructing arm residual forces: mesh points; max error: ', num2str(assertArmTmaxk)])
    end
    if assertMtpTmaxk > 1*10^(-tol_ipopt)
        disp(['Issue when reconstructing mtp residual forces: mesh points; max error: ', num2str(assertMtpTmaxk)])
    end
    
    %% Joint torques and ground reaction forces at collocation points
    Xj_Qs_Qdots_opt             = zeros(d*N,2*nq.all);
    Xj_Qs_Qdots_opt(:,1:2:end)  = q_col_opt_unsc.rad;
    Xj_Qs_Qdots_opt(:,2:2:end)  = qdot_col_opt_unsc.rad;
    Xj_Qdotdots_opt             = qdotdot_col_opt_unsc.rad;
    Foutj_opt = zeros(d*N,nq.all+NGRF+NcalcOrall);
    Tau_passj_opt_all = zeros(d*N,nq.all-nq.abs);
    for i = 1:d*N
        [res] = F1([Xj_Qs_Qdots_opt(i,:)';Xj_Qdotdots_opt(i,:)']);
        Foutj_opt(i,:) = full(res);
        
        Tau_passj_opt.hip.flex.l = full(f_PassiveMoments(k_pass.hip.flex,...
            theta.pass.hip.flex,q_col_opt_unsc.rad(i,jointi.hip_flex.l),...
            qdot_col_opt_unsc.rad(i,jointi.hip_flex.l)));
        Tau_passj_opt.hip.flex.r = full(f_PassiveMoments(k_pass.hip.flex,...
            theta.pass.hip.flex,q_col_opt_unsc.rad(i,jointi.hip_flex.r),...
            qdot_col_opt_unsc.rad(i,jointi.hip_flex.r)));
        Tau_passj_opt.hip.add.l = full(f_PassiveMoments(k_pass.hip.add,...
            theta.pass.hip.add,q_col_opt_unsc.rad(i,jointi.hip_add.l),...
            qdot_col_opt_unsc.rad(i,jointi.hip_add.l)));
        Tau_passj_opt.hip.add.r = full(f_PassiveMoments(k_pass.hip.add,...
            theta.pass.hip.add,q_col_opt_unsc.rad(i,jointi.hip_add.r),...
            qdot_col_opt_unsc.rad(i,jointi.hip_add.r)));
        Tau_passj_opt.hip.rot.l = full(f_PassiveMoments(k_pass.hip.rot,...
            theta.pass.hip.rot,q_col_opt_unsc.rad(i,jointi.hip_rot.l),...
            qdot_col_opt_unsc.rad(i,jointi.hip_rot.l)));
        Tau_passj_opt.hip.rot.r = full(f_PassiveMoments(k_pass.hip.rot,...
            theta.pass.hip.rot,q_col_opt_unsc.rad(i,jointi.hip_rot.r),...
            qdot_col_opt_unsc.rad(i,jointi.hip_rot.r)));
        Tau_passj_opt.knee.l = full(f_PassiveMoments(k_pass.knee,...
            theta.pass.knee,q_col_opt_unsc.rad(i,jointi.knee.l),...
            qdot_col_opt_unsc.rad(i,jointi.knee.l)));
        Tau_passj_opt.knee.r = full(f_PassiveMoments(k_pass.knee,...
            theta.pass.knee,q_col_opt_unsc.rad(i,jointi.knee.r),...
            qdot_col_opt_unsc.rad(i,jointi.knee.r)));
        Tau_passj_opt.ankle.l = full(f_PassiveMoments(k_pass.ankle,...
            theta.pass.ankle,q_col_opt_unsc.rad(i,jointi.ankle.l),...
            qdot_col_opt_unsc.rad(i,jointi.ankle.l)));
        Tau_passj_opt.ankle.r = full(f_PassiveMoments(k_pass.ankle,...
            theta.pass.ankle,q_col_opt_unsc.rad(i,jointi.ankle.r),...
            qdot_col_opt_unsc.rad(i,jointi.ankle.r)));
        Tau_passj_opt.subt.l = full(f_PassiveMoments(k_pass.subt,...
            theta.pass.subt,q_col_opt_unsc.rad(i,jointi.subt.l),...
            qdot_col_opt_unsc.rad(i,jointi.subt.l)));
        Tau_passj_opt.subt.r = full(f_PassiveMoments(k_pass.subt,...
            theta.pass.subt,q_col_opt_unsc.rad(i,jointi.subt.r),...
            qdot_col_opt_unsc.rad(i,jointi.subt.r)));
        Tau_passj_opt.trunk.ext = full(f_PassiveMoments(k_pass.trunk.ext,...
            theta.pass.trunk.ext,q_col_opt_unsc.rad(i,jointi.trunk.ext),...
            qdot_col_opt_unsc.rad(i,jointi.trunk.ext)));
        Tau_passj_opt.trunk.ben = full(f_PassiveMoments(k_pass.trunk.ben,...
            theta.pass.trunk.ben,q_col_opt_unsc.rad(i,jointi.trunk.ben),...
            qdot_col_opt_unsc.rad(i,jointi.trunk.ben)));
        Tau_passj_opt.trunk.rot = full(f_PassiveMoments(k_pass.trunk.rot,...
            theta.pass.trunk.rot,q_col_opt_unsc.rad(i,jointi.trunk.rot),...
            qdot_col_opt_unsc.rad(i,jointi.trunk.rot)));
        
        Tau_passj_opt.sh_flex.l = full(f_passiveTATorques(stiffnessArm, dampingArm, ...
            q_col_opt_unsc.rad(i,jointi.sh_flex.l), qdot_col_opt_unsc.rad(i,jointi.sh_flex.l)));
        Tau_passj_opt.sh_add.l = full(f_passiveTATorques(stiffnessArm, dampingArm, ...
            q_col_opt_unsc.rad(i,jointi.sh_add.l), qdot_col_opt_unsc.rad(i,jointi.sh_add.l)));
        Tau_passj_opt.sh_rot.l = full(f_passiveTATorques(stiffnessArm, dampingArm, ...
            q_col_opt_unsc.rad(i,jointi.sh_rot.l), qdot_col_opt_unsc.rad(i,jointi.sh_rot.l)));
        Tau_passj_opt.sh_flex.r = full(f_passiveTATorques(stiffnessArm, dampingArm, ...
            q_col_opt_unsc.rad(i,jointi.sh_flex.r), qdot_col_opt_unsc.rad(i,jointi.sh_flex.r)));
        Tau_passj_opt.sh_add.r = full(f_passiveTATorques(stiffnessArm, dampingArm, ...
            q_col_opt_unsc.rad(i,jointi.sh_add.r), qdot_col_opt_unsc.rad(i,jointi.sh_add.r)));
        Tau_passj_opt.sh_rot.r = full(f_passiveTATorques(stiffnessArm, dampingArm, ...
            q_col_opt_unsc.rad(i,jointi.sh_rot.r), qdot_col_opt_unsc.rad(i,jointi.sh_rot.r)));
        Tau_passj_opt.elb.l = full(f_passiveTATorques(stiffnessArm, dampingArm, ...
            q_col_opt_unsc.rad(i,jointi.elb.l), qdot_col_opt_unsc.rad(i,jointi.elb.l)));
        Tau_passj_opt.elb.r = full(f_passiveTATorques(stiffnessArm, dampingArm, ...
            q_col_opt_unsc.rad(i,jointi.elb.r), qdot_col_opt_unsc.rad(i,jointi.elb.r)));
        Tau_passj_opt.arm(i,:) = [Tau_passj_opt.sh_flex.l, Tau_passj_opt.sh_add.l, ...
            Tau_passj_opt.sh_rot.l, Tau_passj_opt.sh_flex.r, Tau_passj_opt.sh_add.r, ...
            Tau_passj_opt.sh_rot.r, Tau_passj_opt.elb.l, Tau_passj_opt.elb.r];
        
        Tau_passj_opt.mtp.l = full(f_passiveTATorques(stiffnessMtp, dampingMtp, ...
            q_col_opt_unsc.rad(i,jointi.mtp.l), qdot_col_opt_unsc.rad(i,jointi.mtp.l)));
        Tau_passj_opt.mtp.r = full(f_passiveTATorques(stiffnessMtp, dampingMtp, ...
            q_col_opt_unsc.rad(i,jointi.mtp.r), qdot_col_opt_unsc.rad(i,jointi.mtp.r)));
        Tau_passj_opt.mtp.all(i,:) = [Tau_passj_opt.mtp.l, Tau_passj_opt.mtp.r];
        
        Tau_passj_opt_all(i,:) = [Tau_passj_opt.hip.flex.l,...
            Tau_passj_opt.hip.add.l,Tau_passj_opt.hip.rot.l,...
            Tau_passj_opt.hip.flex.r,Tau_passj_opt.hip.add.r,...
            Tau_passj_opt.hip.rot.r,Tau_passj_opt.knee.l,...
            Tau_passj_opt.knee.r,Tau_passj_opt.ankle.l,...
            Tau_passj_opt.ankle.r,Tau_passj_opt.subt.l,...
            Tau_passj_opt.subt.r,Tau_passj_opt.mtp.all(i,:),...
            Tau_passj_opt.trunk.ext,Tau_passj_opt.trunk.ben,...
            Tau_passj_opt.trunk.rot,Tau_passj_opt.arm(i,:)];
        
    end
    assertArmTmaxj = max(max(abs(Foutj_opt(:,armsi)/scaling.ArmTau - ...
        (a_a_col_opt_unsc + Tau_passj_opt.arm/scaling.ArmTau))));
    assertMtpTmaxj = max(max(abs(Foutj_opt(:,mtpi)/scaling.MtpTau - ...
        (a_mtp_col_opt_unsc + Tau_passj_opt.mtp.all/scaling.MtpTau))));
    if assertArmTmaxj > 1*10^(-tol_ipopt)
        disp(['Issue when reconstructing arm residual forces: collocation points; max error: ', num2str(assertArmTmaxj)])
    end
    if assertMtpTmaxj > 1*10^(-tol_ipopt)
        disp(['Issue when reconstructing mtp residual forces: collocation points; max error: ', num2str(assertMtpTmaxj)])
    end
    
    %% Stride length and width
    % For the stride length we also need the values at the end of the
    % interval so N+1 where states but not controls are defined
    Xk_Qs_Qdots_opt_all = zeros(N+1,2*size(q_opt_unsc_all.rad,2));
    Xk_Qs_Qdots_opt_all(:,1:2:end)  = q_opt_unsc_all.rad;
    Xk_Qs_Qdots_opt_all(:,2:2:end)  = qdot_opt_unsc_all.rad;
    % We just want to extract the positions of the calcaneus origins so we
    % do not really care about Qdotdot that we set to 0
    Xk_Qdotdots_opt_all = zeros(N+1,size(q_opt_unsc_all.rad,2));
    out_res_opt_all = zeros(N+1,nq.all+NGRF+NcalcOrall);
    for i = 1:N+1
        [res] = F1([Xk_Qs_Qdots_opt_all(i,:)';Xk_Qdotdots_opt_all(i,:)']);
        out_res_opt_all(i,:) = full(res);
    end
    % The stride length is the distance covered by the calcaneus origin
    % Right leg
    dist_r = sqrt(f_Jnn3(out_res_opt_all(end,calcOrall.r)-...
        out_res_opt_all(1,calcOrall.r)));
    % Left leg
    dist_l = sqrt(f_Jnn3(out_res_opt_all(end,calcOrall.l)-...
        out_res_opt_all(1,calcOrall.l)));
    % The total stride length is the sum of the right and left stride
    % lengths after a half gait cycle, since we assume symmetry
    StrideLength_opt = full(dist_r + dist_l);
    % The stride width is the medial distance between the calcaneus origins
    StepWidth_opt = full(abs(out_res_opt_all(:,calcOrall.r(3)) - ...
        out_res_opt_all(:,calcOrall.l(3))));
    stride_width_mean = mean(StepWidth_opt);
    stride_width_std = std(StepWidth_opt);
    
    %% Assert average speed
    dist_trav_opt = q_opt_unsc_all.rad(end,jointi.pelvis.tx) - ...
        q_opt_unsc_all.rad(1,jointi.pelvis.tx); % distance traveled
    time_elaps_opt = tf_opt; % time elapsed
    vel_aver_opt = dist_trav_opt/time_elaps_opt;
    % assert_v_tg should be 0
    assert_v_tg = abs(vel_aver_opt-v_tgt);
    if assert_v_tg > 1*10^(-tol_ipopt)
        disp('Issue when reconstructing average speed')
    end
    
    %% Decompose optimal cost
    J_opt           = 0;
    E_cost          = 0;
    A_cost          = 0;
    Arm_cost        = 0;
    Mtp_cost        = 0;
    Qdotdot_cost    = 0;
    Pass_cost       = 0;
    GRF_cost        = 0;
    vA_cost         = 0;
    dFTtilde_cost   = 0;
    QdotdotArm_cost = 0;
    count           = 1;
    h_opt           = tf_opt/N;
    for k=1:N
        for j=1:d
            % Get muscle-tendon lengths, velocities, moment arms
            % Left leg
            qin_l_opt_all = Xj_Qs_Qdots_opt(count,IndexLeft*2-1);
            qdotin_l_opt_all = Xj_Qs_Qdots_opt(count,IndexLeft*2);
            [lMTk_l_opt_all,vMTk_l_opt_all,~] = ...
                f_lMT_vMT_dM(qin_l_opt_all,qdotin_l_opt_all);
            % Right leg
            qin_r_opt_all = Xj_Qs_Qdots_opt(count,IndexRight*2-1);
            qdotin_r_opt_all = Xj_Qs_Qdots_opt(count,IndexRight*2);
            [lMTk_r_opt_all,vMTk_r_opt_all,~] = ...
                f_lMT_vMT_dM(qin_r_opt_all,qdotin_r_opt_all);
            % Both legs
            lMTk_lr_opt_all = ...
                [lMTk_l_opt_all([1:43,47:49],1);lMTk_r_opt_all(1:46,1)];
            vMTk_lr_opt_all = ...
                [vMTk_l_opt_all([1:43,47:49],1);vMTk_r_opt_all(1:46,1)];
            % force equilibirum
            [~,~,Fce_opt_all,Fpass_opt_all,Fiso_opt_all,vMmax_opt_all,...
                massM_opt_all] = f_forceEquilibrium_FtildeState_all_tendon(...
                a_col_opt_unsc(count,:)',FTtilde_col_opt_unsc(count,:)',...
                dFTtilde_col_opt_unsc(count,:)',full(lMTk_lr_opt_all),...
                full(vMTk_lr_opt_all),tensions,aTendon,shift);
            % muscle-tendon kinematics
            [~,lMtilde_opt_all] = f_FiberLength_TendonForce_tendon(...
                FTtilde_col_opt_unsc(count,:)',full(lMTk_lr_opt_all),aTendon,...
                shift);
            [vM_opt_all,~] = f_FiberVelocity_TendonForce_tendon(...
                FTtilde_col_opt_unsc(count,:)',...
                dFTtilde_col_opt_unsc(count,:)',full(lMTk_lr_opt_all),...
                full(vMTk_lr_opt_all),aTendon,shift);
            
            % Bhargava et al. (2004)
            [e_tot_all,~,~,~,~,~] = fgetMetabolicEnergySmooth2004all(...
                a_col_opt_unsc(count,:)',a_col_opt_unsc(count,:)',...
                full(lMtilde_opt_all),...
                full(vM_opt_all),full(Fce_opt_all)',full(Fpass_opt_all)',...
                full(massM_opt_all)',pctsts,full(Fiso_opt_all)',...
                MTparameters_m(1,:)',body_mass,10);            
            e_tot_opt_all = full(e_tot_all)';
            
            % objective function
            J_opt = J_opt + 1/(dist_trav_opt)*(...
                W.E*B(j+1) * (f_J92exp(e_tot_opt_all,exp_E))/body_mass*h_opt + ...
                W.A*B(j+1) * (f_J92(a_col_opt(count,:)))*h_opt +...
                W.ArmE*B(j+1) * (f_J8(e_a_opt(k,:)))*h_opt +...
                W.Mtp*B(j+1) * (f_J2(e_mtp_opt(k,:)))*h_opt +...
                W.Ak*B(j+1) * (f_J23(qdotdot_col_opt(count,residuals_noarmsi)))*h_opt +...
                W.passMom*B(j+1)* (f_J25(Tau_passj_opt_all(count,:)))*h_opt + ...
                W.u*B(j+1) * (f_J92(vA_opt(k,:)))*h_opt + ...
                W.u*B(j+1) * (f_J92(dFTtilde_col_opt(count,:)))*h_opt + ...
                W.u*B(j+1) * (f_J8(qdotdot_col_opt(count,armsi)))*h_opt);
            
            E_cost = E_cost + W.E*B(j+1)*...
                (f_J92exp(e_tot_opt_all,exp_E))/body_mass*h_opt;
            A_cost = A_cost + W.A*B(j+1)*...
                (f_J92(a_col_opt(count,:)))*h_opt;
            Arm_cost = Arm_cost + W.ArmE*B(j+1)*...
                (f_J8(e_a_opt(k,:)))*h_opt;
            Mtp_cost = Mtp_cost + W.Mtp*B(j+1)*...
                (f_J2(e_mtp_opt(k,:)))*h_opt;
            Qdotdot_cost = Qdotdot_cost + W.Ak*B(j+1)*...
                (f_J23(qdotdot_col_opt(count,residuals_noarmsi)))*h_opt;
            Pass_cost = Pass_cost + W.passMom*B(j+1)*...
                (f_J25(Tau_passj_opt_all(count,:)))*h_opt;
            vA_cost = vA_cost + W.u*B(j+1)*...
                (f_J92(vA_opt(k,:)))*h_opt;
            dFTtilde_cost = dFTtilde_cost + W.u*B(j+1)*...
                (f_J92(dFTtilde_col_opt(count,:)))*h_opt;
            QdotdotArm_cost = QdotdotArm_cost + W.u*B(j+1)*...
                (f_J8(qdotdot_col_opt(count,armsi)))*h_opt;
            count = count + 1;
        end
    end
    J_optf = full(J_opt);
    E_costf = full(E_cost);
    A_costf = full(A_cost);
    Arm_costf = full(Arm_cost);
    Mtp_costf = full(Mtp_cost);
    Qdotdot_costf = full(Qdotdot_cost);
    Pass_costf = full(Pass_cost);
    vA_costf = full(vA_cost);
    dFTtilde_costf = full(dFTtilde_cost);
    QdotdotArm_costf = full(QdotdotArm_cost);
    % assertCost should be 0
    assertCost = abs(J_optf - 1/(dist_trav_opt)*(E_costf+A_costf+Arm_costf+...
        Mtp_costf+Qdotdot_costf+Pass_costf+vA_costf+dFTtilde_costf+...
        QdotdotArm_costf));
    assertCost2 = abs(stats.iterations.obj(end) - J_optf);
    if assertCost > 1*10^(-tol_ipopt)
        disp('Issue when reconstructing optimal cost wrt sum of terms')
    end
    if assertCost2 > 1*10^(-tol_ipopt)
        disp('Issue when reconstructing optimal cost wrt stats')
    end
    
    %% Reconstruct full gait cycle
    % We reconstruct the full gait cycle from the simulated half gait cycle
    % Identify heel strike
    threshold = 20; % there is foot-ground contact above the threshold
%     % We adjust some thresholds manually
%     if ww == 17
%         threshold = 31;
%     end
%     
    if exist('HS1','var')
        clear HS1
    end
    % Check if heel strike is on the right side
    phase_tran_tgridi = find(GRFk_opt(:,2)<threshold,1,'last');
    if ~isempty(phase_tran_tgridi)
        if phase_tran_tgridi == N
            temp_idx = find(GRFk_opt(:,2)>threshold,1,'first');
            if ~isempty(temp_idx)
                if temp_idx-1 ~= 0 && ...
                        find(GRFk_opt(temp_idx-1,2)<threshold)
                    phase_tran_tgridi_t = temp_idx;
                    IC1i = phase_tran_tgridi_t;
                    HS1 = 'r';
                end
            else
                IC1i = phase_tran_tgridi + 1;
                HS1 = 'r';
            end
        else
            IC1i = phase_tran_tgridi + 1;
            HS1 = 'r';
        end
    end
    if ~exist('HS1','var')
        % Check if heel strike is on the left side
        phase_tran_tgridi = find(GRFk_opt(:,5)<threshold,1,'last');
        if phase_tran_tgridi == N
            temp_idx = find(GRFk_opt(:,5)>threshold,1,'first');
            if ~isempty(temp_idx)
                if temp_idx-1 ~= 0 && ...
                        find(GRFk_opt(temp_idx-1,5)<threshold)
                    phase_tran_tgridi_t = temp_idx;
                    IC1i = phase_tran_tgridi_t;
                    HS1 = 'l';
                else
                    IC1i = phase_tran_tgridi + 1;
                    HS1 = 'l';
                end
            else
                IC1i = phase_tran_tgridi + 1;
                HS1 = 'l';
            end
        else
            IC1i = phase_tran_tgridi + 1;
            HS1 = 'l';
        end
    end
    
    % GRFk_opt is at mesh points starting from k=2, we thus add 1 to IC1i
    % for the states
    IC1i_c = IC1i;
    IC1i_s = IC1i + 1;
    
%     if isempty(phase_tran_tgridi)
%         disp('No heel strike detected, consider increasing the threshold');
%         continue;
%     end
    
    % Joint kinematics
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
    
    % Qs
    Qs_GC = zeros(N*2,size(q_opt_unsc.deg,2));
    Qs_GC(1:N-IC1i_s+1,:) = q_opt_unsc.deg(IC1i_s:end,:);
    Qs_GC(N-IC1i_s+2:N-IC1i_s+1+N,QsSymA) = q_opt_unsc.deg(1:end,QsSymB);
    Qs_GC(N-IC1i_s+2:N-IC1i_s+1+N,QsOpp) = -q_opt_unsc.deg(1:end,QsOpp);
    Qs_GC(N-IC1i_s+2:N-IC1i_s+1+N,jointi.pelvis.tx) = ...
        q_opt_unsc.deg(1:end,jointi.pelvis.tx) + ...
        q_opt_unsc_all.deg(end,jointi.pelvis.tx);
    Qs_GC(N-IC1i_s+2+N:2*N,:) = q_opt_unsc.deg(1:IC1i_s-1,:);
    Qs_GC(N-IC1i_s+2+N:2*N,jointi.pelvis.tx) = ...
        q_opt_unsc.deg(1:IC1i_s-1,jointi.pelvis.tx) + ...
        2*q_opt_unsc_all.deg(end,jointi.pelvis.tx);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        Qs_GC(:,QsSymA_ptx)  = Qs_GC(:,QsSymB_ptx);
        Qs_GC(:,QsOpp)       = -Qs_GC(:,QsOpp);
    end
    temp_Qs_GC_pelvis_tx = Qs_GC(1,jointi.pelvis.tx);
    Qs_GC(:,jointi.pelvis.tx) = Qs_GC(:,jointi.pelvis.tx)-...
        temp_Qs_GC_pelvis_tx;
    
    % Qdots
    Qdots_GC = zeros(N*2,size(Qs_GC,2));
    Qdots_GC(1:N-IC1i_s+1,:) = qdot_opt_unsc.deg(IC1i_s:end,:);
    Qdots_GC(N-IC1i_s+2:N-IC1i_s+1+N,QsSymA_ptx) = ...
        qdot_opt_unsc.deg(1:end,QsSymB_ptx);
    Qdots_GC(N-IC1i_s+2:N-IC1i_s+1+N,QsOpp) = ...
        -qdot_opt_unsc.deg(1:end,QsOpp);
    Qdots_GC(N-IC1i_s+2+N:2*N,:) = qdot_opt_unsc.deg(1:IC1i_s-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        Qdots_GC(:,QsSymA_ptx) = Qdots_GC(:,QsSymB_ptx);
        Qdots_GC(:,QsOpp) = -Qdots_GC(:,QsOpp);
    end
    
    % Qdotdots
    Qdotdots_GC = zeros(N*2,size(Qs_opt,2));
    Qdotdots_GC(1:N-IC1i_c+1,:) = Xk_Qdotdots_opt(IC1i_c:end,:);
    Qdotdots_GC(N-IC1i_c+2:N-IC1i_c+1+N,QsSymA_ptx) = ...
        Xk_Qdotdots_opt(1:end,QsSymB_ptx);
    Qdotdots_GC(N-IC1i_c+2:N-IC1i_c+1+N,QsOpp) = ...
        -Xk_Qdotdots_opt(1:end,QsOpp);
    Qdotdots_GC(N-IC1i_c+2+N:2*N,:) = Xk_Qdotdots_opt(1:IC1i_c-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        Qdotdots_GC(:,QsSymA_ptx) = Qdotdots_GC(:,QsSymB_ptx);
        Qdotdots_GC(:,QsOpp) = -Qdotdots_GC(:,QsOpp);
    end
    
    % Ground reaction forces
    GRFs_opt = zeros(N*2,NGRF);
    GRFs_opt(1:N-IC1i_c+1,:) = GRFk_opt(IC1i_c:end,1:6);
    GRFs_opt(N-IC1i_c+2:N-IC1i_c+1+N,:) = GRFk_opt(1:end,[4:6,1:3]);
    GRFs_opt(N-IC1i_c+2:N-IC1i_c+1+N,[3,6]) = ...
        -GRFs_opt(N-IC1i_c+2:N-IC1i_c+1+N,[3,6]);
    GRFs_opt(N-IC1i_c+2+N:2*N,:) = GRFk_opt(1:IC1i_c-1,1:6);
    GRFs_opt = GRFs_opt./(body_weight/100);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        GRFs_opt(:,[4:6,1:3]) = GRFs_opt(:,:);
        GRFs_opt(:,[3,6]) = -GRFs_opt(:,[3,6]);
    end
    
    % Joint torques
    Ts_opt = zeros(N*2,size(Qs_opt,2));
    Ts_opt(1:N-IC1i_c+1,1:nq.all) = Foutk_opt(IC1i_c:end,1:nq.all);
    Ts_opt(N-IC1i_c+2:N-IC1i_c+1+N,QsSymA_ptx) = Foutk_opt(1:end,QsSymB_ptx);
    Ts_opt(N-IC1i_c+2:N-IC1i_c+1+N,QsOpp) = -Foutk_opt(1:end,QsOpp);
    Ts_opt(N-IC1i_c+2+N:2*N,1:nq.all) = Foutk_opt(1:IC1i_c-1,1:nq.all);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        Ts_opt(:,QsSymA_ptx) = Ts_opt(:,QsSymB_ptx);
        Ts_opt(:,QsOpp) = -Ts_opt(:,QsOpp);
    end
    Ts_opt = Ts_opt./body_mass;
    
    % Muscle-Tendon Forces
    orderMusInv = [NMuscle/2+1:NMuscle,1:NMuscle/2];
    FTtilde_GC = zeros(N*2,NMuscle);
    FTtilde_GC(1:N-IC1i_s+1,:) = FTtilde_opt_unsc(IC1i_s:end,:);
    FTtilde_GC(N-IC1i_s+2:N-IC1i_s+1+N,:) = ...
        FTtilde_opt_unsc(1:end,orderMusInv);
    FTtilde_GC(N-IC1i_s+2+N:2*N,:) = FTtilde_opt_unsc(1:IC1i_s-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        FTtilde_GC(:,:) = FTtilde_GC(:,orderMusInv);
    end
    
    % Muscle activations
    Acts_GC = zeros(N*2,NMuscle);
    Acts_GC(1:N-IC1i_s+1,:) = a_opt_unsc(IC1i_s:end,:);
    Acts_GC(N-IC1i_s+2:N-IC1i_s+1+N,:) = a_opt_unsc(1:end,orderMusInv);
    Acts_GC(N-IC1i_s+2+N:2*N,:) = a_opt_unsc(1:IC1i_s-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        Acts_GC(:,:) = Acts_GC(:,orderMusInv);
    end
    
    % Time derivative of muscle-tendon force
    dFTtilde_GC = zeros(N*2,NMuscle);
    dFTtilde_GC(1:N-IC1i_c+1,:) = dFTtilde_opt_unsc(IC1i_c:end,:);
    dFTtilde_GC(N-IC1i_c+2:N-IC1i_c+1+N,:) = ...
        dFTtilde_opt_unsc(1:end,orderMusInv);
    dFTtilde_GC(N-IC1i_c+2+N:2*N,:) = dFTtilde_opt_unsc(1:IC1i_c-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        dFTtilde_GC(:,:) = dFTtilde_GC(:,orderMusInv);
    end
    
    % Muscle excitations
    vA_GC = zeros(N*2,NMuscle);
    vA_GC(1:N-IC1i_c+1,:) = vA_opt_unsc(IC1i_c:end,:);
    vA_GC(N-IC1i_c+2:N-IC1i_c+1+N,:) = vA_opt_unsc(1:end,orderMusInv);
    vA_GC(N-IC1i_c+2+N:2*N,:) = vA_opt_unsc(1:IC1i_c-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        vA_GC(:,:) = vA_GC(:,orderMusInv);
    end
    e_GC = computeExcitationRaasch(Acts_GC,vA_GC,...
        ones(1,NMuscle)*tdeact,ones(1,NMuscle)*tact);
    
    % Arm activations
    orderArmInv = [jointi.sh_flex.r:jointi.sh_rot.r,...
        jointi.sh_flex.l:jointi.sh_rot.l,...
        jointi.elb.r,jointi.elb.l]-jointi.sh_flex.l+1;
    a_a_GC = zeros(N*2,nq.arms);
    a_a_GC(1:N-IC1i_s+1,:) = a_a_opt_unsc(IC1i_s:end,:);
    a_a_GC(N-IC1i_s+2:N-IC1i_s+1+N,:) = a_a_opt_unsc(1:end,orderArmInv);
    a_a_GC(N-IC1i_s+2+N:2*N,:) = a_a_opt_unsc(1:IC1i_s-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        a_a_GC(:,:) = a_a_GC(:,orderArmInv);
    end
    
    % Mtp activations
    orderMtpInv = [jointi.mtp.r,jointi.mtp.l]-jointi.mtp.l+1;
    a_mtp_GC = zeros(N*2,nq.mtp);
    a_mtp_GC(1:N-IC1i_s+1,:) = a_mtp_opt_unsc(IC1i_s:end,:);
    a_mtp_GC(N-IC1i_s+2:N-IC1i_s+1+N,:) = a_mtp_opt_unsc(1:end,orderMtpInv);
    a_mtp_GC(N-IC1i_s+2+N:2*N,:) = a_mtp_opt_unsc(1:IC1i_s-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        a_mtp_GC(:,:) = a_mtp_GC(:,orderMtpInv);
    end
    
    % Arm excitations
    e_a_GC = zeros(N*2,nq.arms);
    e_a_GC(1:N-IC1i_c+1,:) = e_a_opt_unsc(IC1i_c:end,:);
    e_a_GC(N-IC1i_c+2:N-IC1i_c+1+N,:) = e_a_opt_unsc(1:end,orderArmInv);
    e_a_GC(N-IC1i_c+2+N:2*N,:) = e_a_opt_unsc(1:IC1i_c-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        e_a_GC(:,:) = e_a_GC(:,orderArmInv);
    end
    
    % Mtp excitations
    e_mtp_GC = zeros(N*2,nq.mtp);
    e_mtp_GC(1:N-IC1i_c+1,:) = e_mtp_opt_unsc(IC1i_c:end,:);
    e_mtp_GC(N-IC1i_c+2:N-IC1i_c+1+N,:) = e_mtp_opt_unsc(1:end,orderMtpInv);
    e_mtp_GC(N-IC1i_c+2+N:2*N,:) = e_mtp_opt_unsc(1:IC1i_c-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        e_mtp_GC(:,:) = e_mtp_GC(:,orderMtpInv);
    end
    
    % ExoTorques
    T_exo_GC = zeros(N*2,2);
    T_exo_GC(1:N-IC1i_c+1,:) = ExoVect([1 2],IC1i_c:end)';
    T_exo_GC(N-IC1i_c+2:N-IC1i_c+1+N,:) = ExoVect([2 1],1:end)';
    T_exo_GC(N-IC1i_c+2+N:2*N,:) = ExoVect([1 2],1:IC1i_c-1)';    
    
    % Passive joint torques
    Tau_pass_opt_inv = [jointi.hip_flex.r:jointi.hip_rot.r,...
        jointi.hip_flex.l:jointi.hip_rot.l,...
        jointi.knee.r,jointi.knee.l,jointi.ankle.r,jointi.ankle.l,...
        jointi.subt.r,jointi.subt.l,jointi.mtp.r,jointi.mtp.l,...
        jointi.trunk.ext:jointi.trunk.rot,...
        jointi.sh_flex.r:jointi.sh_rot.r,...
        jointi.sh_flex.l:jointi.sh_rot.l,...
        jointi.elb.r,jointi.elb.l]-jointi.hip_flex.l+1;
    Tau_pass_opt_GC = zeros(N*2,nq.all-nq.abs);
    Tau_pass_opt_GC(1:N-IC1i_c+1,:) = Tau_passk_opt_all(IC1i_c:end,:);
    Tau_pass_opt_GC(N-IC1i_c+2:N-IC1i_c+1+N,:) = ...
        Tau_passk_opt_all(1:end,Tau_pass_opt_inv);
    Tau_pass_opt_GC(N-IC1i_c+2+N:2*N,:) = Tau_passk_opt_all(1:IC1i_c-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        Tau_pass_opt_GC(:,Tau_pass_opt_inv) = Tau_pass_opt_GC(:,:);
    end
    
    % Create .mot file for OpenSim GUI
    q_opt_GUI_GC = zeros(2*N,1+nq.all+2);
    q_opt_GUI_GC(1:N-IC1i_s+1,1) = tgrid(:,IC1i_s:end-1)';
    q_opt_GUI_GC(N-IC1i_s+2:N-IC1i_s+1+N,1)  = tgrid(:,1:end-1)' + tgrid(end);
    q_opt_GUI_GC(N-IC1i_s+2+N:2*N,1) = tgrid(:,1:IC1i_s-1)' + 2*tgrid(end);
    q_opt_GUI_GC(:,2:end-2) = Qs_GC;
    q_opt_GUI_GC(:,end-1:end) = 1.51*180/pi*ones(2*N,2); % pro_sup (locked)
    q_opt_GUI_GC(:,1) = q_opt_GUI_GC(:,1)-q_opt_GUI_GC(1,1);
    if writeIKmotion
        pathOpenSim = [pathRepo,'/OpenSim'];
        addpath(genpath(pathOpenSim));
        JointAngle.labels = {'time','pelvis_tilt','pelvis_list',...
            'pelvis_rotation','pelvis_tx','pelvis_ty','pelvis_tz',...
            'hip_flexion_l','hip_adduction_l','hip_rotation_l',...
            'hip_flexion_r','hip_adduction_r','hip_rotation_r',...
            'knee_angle_l','knee_angle_r','ankle_angle_l','ankle_angle_r',...
            'subtalar_angle_l','subtalar_angle_r','mtp_angle_l','mtp_angle_r',...
            'lumbar_extension','lumbar_bending','lumbar_rotation',...
            'arm_flex_l','arm_add_l','arm_rot_l',...
            'arm_flex_r','arm_add_r','arm_rot_r',...
            'elbow_flex_l','elbow_flex_r',...
            'pro_sup_l','pro_sup_r'};
        % Two gait cycles
        % Joint angles
        q_opt_GUI_GC_2 = [q_opt_GUI_GC;q_opt_GUI_GC];
        q_opt_GUI_GC_2(2*N+1:4*N,1) = q_opt_GUI_GC_2(2*N+1:4*N,1) + ...
            q_opt_GUI_GC_2(end,1) + ...
            q_opt_GUI_GC_2(end,1)-q_opt_GUI_GC_2(end-1,1);
        q_opt_GUI_GC_2(2*N+1:4*N,jointi.pelvis.tx+1) = ...
            q_opt_GUI_GC_2(2*N+1:4*N,jointi.pelvis.tx+1) + ...
            2*q_opt_unsc_all.deg(end,jointi.pelvis.tx);
        % Muscle activations (to have muscles turning red when activated).
        Acts_GC_GUI = [Acts_GC;Acts_GC];
        % Combine data joint angles and muscle activations
        JointAngleMuscleAct.data = [q_opt_GUI_GC_2,Acts_GC_GUI];
        % Get muscle labels
        muscleNamesAll = cell(1,NMuscle);
        for i = 1:NMuscle/2
            muscleNamesAll{i} = [muscleNames{i}(1:end-2),'_l'];
            muscleNamesAll{i+NMuscle/2} = [muscleNames{i}(1:end-2),'_r'];
        end
        % Combine labels joint angles and muscle activations
        JointAngleMuscleAct.labels = JointAngle.labels;
        for i = 1:NMuscle
            JointAngleMuscleAct.labels{i+size(q_opt_GUI_GC_2,2)} = ...
                [muscleNamesAll{i},'/activation'];
        end
        OutFolder = fullfile(pathRepo,'Results',S.ResultsFolder);
        filenameJointAngles = fullfile(OutFolder,[S.savename '.mot']); 
        write_motionFile(JointAngleMuscleAct, filenameJointAngles);
    end
    
    %% Metabolic cost of transport for a gait cycle
    Qs_opt_rad = Qs_GC;
    Qs_opt_rad(:,roti) = Qs_opt_rad(:,roti).*pi/180;
    qdot_opt_GC_rad = Qdots_GC;
    qdot_opt_GC_rad(:,roti)= qdot_opt_GC_rad(:,roti).*pi/180;
    % Pre-allocations
    e_mo_opt = zeros(2*N,1);
    e_mo_optb = zeros(2*N,1);
    vMtilde_opt_all = zeros(2*N, NMuscle);
    lMtilde_opt_all = zeros(2*N, NMuscle);
    metab_Etot = zeros(2*N, NMuscle);
    metab_Adot = zeros(2*N, NMuscle);
    metab_Mdot = zeros(2*N, NMuscle);
    metab_Sdot = zeros(2*N, NMuscle);
    metab_Wdot = zeros(2*N, NMuscle);
    FT_opt     = zeros(2*N, NMuscle);
    for nn = 1:2*N
        % Get muscle-tendon lengths, velocities, moment arms
        % Left leg
        qin_l_opt = Qs_opt_rad(nn,IndexLeft);
        qdotin_l_opt = qdot_opt_GC_rad(nn,IndexLeft);
        [lMTk_l_opt,vMTk_l_opt,~] = f_lMT_vMT_dM(qin_l_opt,qdotin_l_opt);
        % Right leg
        qin_r_opt = Qs_opt_rad(nn,IndexRight);
        qdotin_r_opt = qdot_opt_GC_rad(nn,IndexRight);
        [lMTk_r_opt,vMTk_r_opt,~] = f_lMT_vMT_dM(qin_r_opt,qdotin_r_opt);        
        % Both legs
        lMTk_lr_opt     = [lMTk_l_opt([1:43,47:49],1);lMTk_r_opt(1:46,1)];
        vMTk_lr_opt     = [vMTk_l_opt([1:43,47:49],1);vMTk_r_opt(1:46,1)];        
        % force equilibrium
        [Hilldiff_optt,FT_optt,Fce_optt,Fpass_optt,Fiso_optt,...
            vMmax_optt,massM_optt] = f_forceEquilibrium_FtildeState_all_tendon(...
            Acts_GC(nn,:)',FTtilde_GC(nn,:)',dFTtilde_GC(nn,:)',full(lMTk_lr_opt),...
            full(vMTk_lr_opt),tensions,aTendon,shift);
        % fiber kinematics
        [~,lMtilde_opt] = f_FiberLength_TendonForce_tendon(...
            FTtilde_GC(nn,:)',full(lMTk_lr_opt),aTendon,shift);
        lMtilde_opt_all(nn,:) = full(lMtilde_opt)';
        [vM_opt,vMtilde_opt] = f_FiberVelocity_TendonForce_tendon(FTtilde_GC(nn,:)',...
            dFTtilde_GC(nn,:)',full(lMTk_lr_opt),full(vMTk_lr_opt),...
            aTendon,shift);
        vMtilde_opt_all(nn,:) = full(vMtilde_opt)';
        % Bhargava et al. (2004)
        [energy_total,Adot,Mdot,Sdot,Wdot,e_mot] = ...
            fgetMetabolicEnergySmooth2004all(Acts_GC(nn,:)',...
            Acts_GC(nn,:)',full(lMtilde_opt),full(vM_opt),...
            full(Fce_optt),full(Fpass_optt),full(massM_optt),pctsts,...
            full(Fiso_optt),MTparameters_m(1,:)',body_mass,10);
        e_mo_opt(nn,:) = full(e_mot)';
        e_mo_optb(nn,:) = full(e_mot)';
        metab_Etot(nn,:) = full(energy_total)';
        metab_Adot(nn,:) = full(Adot)';
        metab_Mdot(nn,:) = full(Mdot)';
        metab_Sdot(nn,:) = full(Sdot)';
        metab_Wdot(nn,:) = full(Wdot)';
        FT_opt(nn,:)     = full(FT_optt)';
    end
    % Get COT
    dist_trav_opt_GC = Qs_opt_rad(end,jointi.pelvis.tx) - ...
        Qs_opt_rad(1,jointi.pelvis.tx); % distance traveled
    time_GC = q_opt_GUI_GC(:,1);
    e_mo_opt_tr = trapz(time_GC,e_mo_opt);
    e_mo_opt_trb = trapz(time_GC,e_mo_optb);
    % Cost of transport: J/kg/m
    % Energy model used for optimization
    COT_GC_InOpt = e_mo_opt_tr/body_mass/dist_trav_opt_GC;
    % Energy model from Bhargava et al. (2004)
    COT_GC = e_mo_opt_trb/body_mass/dist_trav_opt_GC;
    
    %% Save results    
   
    % To Do:
    % 1) export casadi function to speed-up?
    % 2) 
    
    % Structure Results_all
    R.t_step    = tgrid;
    R.tf_step   = tgrid(end);
    R.t         = q_opt_GUI_GC(:,1);
    R.tend      = q_opt_GUI_GC(end,1) - q_opt_GUI_GC(1,1);
    R.Qs        = Qs_GC;
    R.Qdots     = Qdots_GC;
    R.GRFs      = GRFs_opt;
    R.Ts        = Ts_opt;
    R.Tid       = Ts_opt.*body_mass;
    R.a         = Acts_GC;
    R.e         = e_GC;
    R.COT       = COT_GC;
    R.StrideLength = StrideLength_opt;
    R.StepWidth = stride_width_mean;
    R.vMtilde   = vMtilde_opt_all;
    R.lMtilde   = lMtilde_opt_all;
    R.MetabB.Etot = metab_Etot;
    R.MetabB.Adot = metab_Adot;
    R.MetabB.Mdot = metab_Mdot;
    R.MetabB.Sdot = metab_Sdot;
    R.MetabB.Wdot = metab_Wdot;
    R.ExoControl  = ExoControl;
    R.S           = S;  % settings for post processing
    R.Sopt        = Sopt; % original settings used to solve the OCP
    R.body_mass   = body_mass;
    R.a_arm       = a_a_GC;
    R.e_arm       = e_a_GC;
    R.a_mtp       = a_mtp_GC;
    R.e_mtp       = e_mtp_GC;
    R.FT          = FT_opt;
    R.TPass       = Tau_pass_opt_GC;
    R.dt          = nanmean(diff(R.t));
    R.T_exo       = T_exo_GC;
    R.dt_exoShift = IC1i_c.*R.dt;
    
    % header information
    R.colheaders.joints = joints;
    R.colheaders.GRF = {'fore_aft_r','vertical_r',...
        'lateral_r','fore_aft_l','vertical_l','lateral_l'};
    for i = 1:NMuscle/2
        R.colheaders.muscles{i} = ...
            [muscleNames{i}(1:end-2),'_l'];
        R.colheaders.muscles{i+NMuscle/2} = ...
            [muscleNames{i}(1:end-2),'_r'];
    end

    % script information
    R.info.script = 'f_PredSim_PoggenSee2020.m';

    % Save data
    OutFolder = fullfile(pathRepo,'Results',S.ResultsFolder);
    FilenameAnalysis = fullfile(OutFolder,[S.savename '_pp.mat']); 
    save(FilenameAnalysis,'R');
end
