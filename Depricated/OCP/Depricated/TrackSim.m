%%  Three-dimensional muscle-driven tracking simulations of human gaits
%   with optimization of contact model parameters
%
% Author: Antoine Falisse
% Date: 1/7/2019
%
% DEPENDENCY: please install CasADi (https://web.casadi.org/)
%
clear all;
clc
close all;

%% User inputs
% This script can be run to solve the optimal control problems but also to
% analyze and process the results. The user does not need to re-run the
% optimal control problems to analyze the results. Therefore, the user can
% select some settings beforehand through the variable num_set. For
% example, when num_set(1) is set to 0, the script will not run the
% optimal control problem. Here is a brief description of num_set:
% num_set(1): set to 1 to solve problem
% num_set(2): set to 1 to analyze results
% num_set(3): set to 1 to load results
% num_set(4): set to 1 to save results
% num_set(5): set to 1 to visualize guess-bounds 
% num_set(6): set to 1 to write .mot file
% Note that you should re-run the simulations to write out the .mot files
% and visualize the results in the OpenSim GUI.

num_set = [1,1,0,1,0,1]; % This configuration solves the problem
% num_set = [0,1,1,0,0,1]; % This configuration analyzes the results

% The variable settings in the following section will set some parameters 
% of the optimization problem. Through the variable idx_ww, the user can 
% select which row of parameters will be used.
idx_ww = 1; % Index row in matrix settings

%% Settings
import casadi.*
subject = 'subject1';

solveProblem    = num_set(1); % set to 1 to solve problem
analyseResults  = num_set(2); % set to 1 to analyze results
loadResults     = num_set(3); % set to 1 to load results
saveResults     = num_set(4); % set to 1 to save sens. results
checkBoundsIG   = num_set(5); % set to 1 to visualize guess-bounds 
writeIKmotion   = num_set(6); % set to 1 to write .mot file

% settings describes the parameters used in the optimal control problem.
% settings(1): weight joint kinematics
% settings(2): weight ground reaction forces
% settings(3): weight ground reaction moments
% settings(4): weight joint kinetics
% settings(5): weight muscle activations
% settings(6): allowed deviation from generic locations
% settings(7): allowed deviation from generic radii
% settings(8): number of mesh intervals
% settings(9): NLP error tolerance: 1*10^(-settings(9)).
settings = [10,1,1,10,1,25,50,50,4];

%% Select settings
for www = 1:length(idx_ww)
ww = idx_ww(www);
% Variable parameters
W.Qs        = settings(ww,1);  % weight joint kinematics
W.GRF       = settings(ww,2);  % weight ground reaction forces
W.GRM       = settings(ww,3);  % weight ground reaction moments
W.ID_act    = settings(ww,4);  % weight joint kinetics
W.a         = settings(ww,5);  % weight muscle activations
dev_cm.loc  = settings(ww,6);  % allowed deviation from generic locations
dev_cm.rad  = settings(ww,7);  % allowed deviation from generic radii
N           = settings(ww,8);  % number of mesh intervals
tol_ipopt   = settings(ww,9);  % NLP error tolerance: 1*10^(-settings(9))
% Fixed parameter
W.u = 0.001;
% Identifiers for experimental data
nametrial.id    = 'gait_1'; % Experimental walking trial to track
nametrial.ID    = ['ID_',nametrial.id];
nametrial.GRF   = ['GRF_',nametrial.id];
nametrial.IK    = ['IK_',nametrial.id];
switch nametrial.id
    case 'gait_1'
        time_opt = [3.45,4.25];
end  
% The filename used to save the results depends on the settings 
savename = ['_c',num2str(ww),'_id_',nametrial.id,'_N',num2str(N),...
    '_T',num2str(tol_ipopt),'_Q',num2str(W.Qs),'_F',num2str(W.GRF),...
    '_M',num2str(W.GRM),'_A',num2str(W.a),'_devL',num2str(dev_cm.loc),...
    '_devR',num2str(dev_cm.rad),'_ID',num2str(W.ID_act)];

%% Load external functions
% The external function performs inverse dynamics through the
% OpenSim/Simbody C++ API. This external function is compiled as a dll from
% which we create a Function instance using CasADi in MATLAB.
pathmain = pwd;
%   
% We use different external functions. A first external function extracts 
% several parameters of the bodies to which the contact spheres are attached.
% The contact forces are then computed in MATLAB and are inputs of the
% second external function in which the skeleton dynamics is described. The
% motivation for this decoupling is to limit the number of times we need to
% build the model. By defining the contact model in MATLAB, we only need to
% build the model once per external function whereas keeping the contact
% model in the external function would require re-building the model during
% the optimization.
[pathRepo,~,~] = fileparts(pathmain);
pathExternalFunctions = [pathRepo,'/ExternalFunctions'];
% Loading external functions. 
setup.derivatives =  'AD'; % Algorithmic differentiation
if ispc 
    switch setup.derivatives
        case 'AD'     
            cd(pathExternalFunctions);
            F1 = external('F','TrackSim_1.dll'); 
            F2 = external('F','TrackSim_2.dll'); 
    end
elseif ismac
    switch setup.derivatives
        case 'AD'     
            cd(pathExternalFunctions);
            F1 = external('F','TrackSim_1.dylib'); 
            F2 = external('F','TrackSim_2.dylib'); 
    end
else
    disp('Platform not supported')
end   
cd(pathmain);
% This is an example of how to call an external function with some
% numerical values.
% vec1 = -ones(141,1);
% res1 = full(F2(vec1));

%% Indices external function
% External function: F1
calcn.l.omega = 1:3;
calcn.l.v_lin = 4:6;
calcn.r.omega = 7:9;
calcn.r.v_lin = 10:12;
calcn.l.pos = 13:15;
calcn.r.pos = 16:18;
calcn.TR.R.l = 19:27;
calcn.TR.T.l = 28:30;
calcn.TR.R.r = 31:39;
calcn.TR.T.r = 40:42;
toes.l.omega = 43:45;
toes.l.v_lin = 46:48;
toes.r.omega = 49:51;
toes.r.v_lin = 52:54;
toes.l.pos = 55:57;
toes.r.pos = 58:60;
toes.TR.R.l = 61:69;
toes.TR.T.l = 70:72;
toes.TR.R.r = 73:81;
toes.TR.T.r = 82:84;
% External function: F2
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
jointi.trunk.ext    = 19;
jointi.trunk.ben    = 20;
jointi.trunk.rot    = 21;
jointi.sh_flex.l    = 22;
jointi.sh_add.l     = 23;
jointi.sh_rot.l     = 24;
jointi.sh_flex.r    = 25;
jointi.sh_add.r     = 26;
jointi.sh_rot.r     = 27;
jointi.elb.l        = 28;
jointi.elb.r        = 29;
% Vectors of indices for later use
residualsi          = jointi.pelvis.tilt:jointi.elb.r; % all 
ground_pelvisi      = jointi.pelvis.tilt:jointi.pelvis.tz; % ground-pelvis
trunki              = jointi.trunk.ext:jointi.trunk.rot; % trunk
armsi               = jointi.sh_flex.l:jointi.elb.r; % arms
residuals_acti      = jointi.hip_flex.l:jointi.elb.r; % all but gr-pelvis
residual_bptyi      = [jointi.pelvis.tilt:jointi.pelvis.tx,...
    jointi.pelvis.tz:jointi.elb.r]; % all but pelvis_ty
% Number of degrees of freedom for later use
nq.all              = length(residualsi); % all 
nq.abs              = length(ground_pelvisi); % ground-pelvis
nq.act              = nq.all-nq.abs;% all but ground-pelvis
nq.trunk            = length(trunki); % trunk
nq.arms             = length(armsi); % arms
nq.leg              = 9; % #joints needed for polynomials
Qsi                 = 1:2:2*nq.all; % indices Qs only
% Second, GRFs
GRFi.r              = 30:32;
GRFi.l              = 33:35;
GRFi.all            = [GRFi.r,GRFi.l];
nGRF                = length(GRFi.all);
% Third, GRMs
GRMi.r              = 36:38;
GRMi.l              = 39:41;
GRMi.all            = [GRMi.r,GRMi.l];
nGRM                = length(GRMi.all);
% Number contact model parameters
np = 18;  

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
load([pathmusclemodel,'/MTparameters_',subject,'.mat']);
MTparameters_m = [MTparameters(:,musi),MTparameters(:,musi)];
% Indices of the muscles actuating the different joints for later use
pathpolynomial = [pathRepo,'/Polynomials'];
addpath(genpath(pathpolynomial));
tl = load([pathpolynomial,'/muscle_spanning_joint_INFO_',subject,'.mat']);
[~,mai] = MomentArmIndices(muscleNames(1:end-3),...
    tl.muscle_spanning_joint_INFO(1:end-3,:));

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

%% Contact model parameters
% We optimize the locations of the contact spheres (x and z coordinates) 
% and the radii. The other parameters are fixed.
% Indices variable parameters
loci.s1.r.x = 1;
loci.s1.r.z = 2;
loci.s2.r.x = 3;
loci.s2.r.z = 4;
loci.s3.r.x = 5;
loci.s3.r.z = 6;
loci.s4.r.x = 7;
loci.s4.r.z = 8;
loci.s5.r.x = 9;
loci.s5.r.z = 10;
loci.s6.r.x = 11;
loci.s6.r.z = 12;
radi.s1 = 13;
radi.s2 = 14;
radi.s3 = 15;
radi.s4 = 16;
radi.s5 = 17;
radi.s6 = 18;
% Fixed parameters
dissipation = 2;
normal = [0,1,0];
transitionVelo = 0.2;
staticFriction = 0.8;
dynamicFriction = 0.8;
viscousFriction = 0.5; 
stif = 1000000;  
loc.s1.y    = -0.021859;
loc.s2.y    = -0.021859;
loc.s3.y    = -0.021859;
loc.s4.y    = -0.0214476;
loc.s5.y    = -0.021859;
loc.s6.y    = -0.0214476;

%% CasADi functions
% We create several CasADi functions for later use
pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
addpath(genpath(pathCasADiFunctions));
pathContactModel = [pathRepo,'/Contact'];
addpath(genpath(pathContactModel));
% We load some variables for the polynomial approximations
load([pathpolynomial,'/muscle_spanning_joint_INFO_',subject,'.mat']);
load([pathpolynomial,'/MuscleInfo_',subject,'.mat']);
% For the polynomials, we want all independent muscles. So we do not need 
% the muscles from both legs, since we assume bilateral symmetry, but want
% all muscles from the back (indices 47:49).
musi_pol = [musi,47,48,49];
NMuscle_pol = NMuscle/2+3;
CasADiFunctions_tracking

%% Passive joint torques
% We extract the parameters for the passive torques of the lower limbs and
% the trunk
pathPassiveMoments = [pathRepo,'/PassiveMoments'];
addpath(genpath(pathPassiveMoments));
PassiveMomentsData

%% Experimental data
pathData = [pathRepo,'/OpenSimModel/',subject];
joints = {'pelvis_tilt','pelvis_list','pelvis_rotation','pelvis_tx',...
    'pelvis_ty','pelvis_tz','hip_flexion_l','hip_adduction_l',...
    'hip_rotation_l','hip_flexion_r','hip_adduction_r','hip_rotation_r',...
    'knee_angle_l','knee_angle_r','ankle_angle_l','ankle_angle_r',...
    'subtalar_angle_l','subtalar_angle_r',...
    'lumbar_extension','lumbar_bending','lumbar_rotation',...
    'arm_flex_l','arm_add_l','arm_rot_l',...
    'arm_flex_r','arm_add_r','arm_rot_r',...
    'elbow_flex_l','elbow_flex_r'};
pathVariousFunctions = [pathRepo,'/VariousFunctions'];
addpath(genpath(pathVariousFunctions));
% Extract joint kinematics
pathIK = [pathData,'/IK/',nametrial.IK,'.mat'];
Qs = getIK(pathIK,joints);
% Extract ground reaction forces and moments
pathGRF = [pathData,'/GRF/',nametrial.GRF,'.mat'];
GRF = getGRF(pathGRF);
% Extract joint kinetics
pathID = [pathData,'/ID/',nametrial.ID,'.mat'];
ID = getID(pathID,joints);
% Interpolation experimental data
time_expi.ID(1) = find(round(ID.time,4) == time_opt(1));
time_expi.ID(2) = find(round(ID.time,4) == time_opt(2));
time_expi.GRF(1) = find(round(GRF.time,4) == time_opt(1));
time_expi.GRF(2) = find(round(GRF.time,4) == time_opt(2));
step = (ID.time(time_expi.ID(2))-ID.time(time_expi.ID(1)))/(N-1);
interval = time_opt(1):step:time_opt(2);
ID.allinterp = interp1(ID.all(:,1),ID.all,interval);
Qs.allinterpfilt = interp1(Qs.allfilt(:,1),Qs.allfilt,interval);
GRF.val.allinterp = interp1(round(GRF.val.all(:,1),4),...
    GRF.val.all,round(interval,4));
GRF.MorGF.allinterp = interp1(round(GRF.MorGF.all(:,1),4),...
    GRF.MorGF.all,round(interval,4));
GRF.pos.allinterp = interp1(round(GRF.pos.all(:,1),4),...
    GRF.pos.all,round(interval,4));
GRF.Mcop.allinterp = interp1(round(GRF.Mcop.all(:,1),4),...
    GRF.Mcop.all,round(interval,4));

%% Bounds
pathBounds = [pathRepo,'/Bounds'];
addpath(genpath(pathBounds));
[bounds,scaling] = getBounds_tracking(Qs,NMuscle,nq,jointi,dev_cm,GRF);

%% Initial guess
pathIG = [pathRepo,'/IG'];
addpath(genpath(pathIG));
% Data-informed initial guess
guess = getGuess_DI_tracking(Qs,nq,N,NMuscle,jointi,scaling);
% This allows visualizing the initial guess and the bounds
if checkBoundsIG
    pathPlots = [pathRepo,'/Plots'];
    addpath(genpath(pathPlots));
    plot_BoundsVSInitialGuess_tracking
end

%% Formulate the NLP
if solveProblem
    % Start with an empty NLP
    w   = {}; % design variables
    w0  = []; % initial guess for design variables
    lbw = []; % lower bounds for design variables
    ubw = []; % upper bounds for design variables
    J   = 0;  % initial value of cost function
    g   = {}; % constraints
    lbg = []; % lower bounds for constraints
    ubg = []; % upper bounds for constraints
    % Define static parameters
    % Contact model parameters   
    paramsCM        = MX.sym('paramsCM',np);
    w               = [w {paramsCM}];
    lbw             = [lbw; bounds.params.lower'];
    ubw             = [ubw; bounds.params.upper'];
    w0              = [w0;  guess.params'];      
    paramsCM_nsc    = f_nsc18(paramsCM,scaling.params.v,scaling.params.r); 
    % Define states at first mesh point
    % Muscle activations
    a0              = MX.sym('a0',NMuscle);
    w               = [w {a0}];
    lbw             = [lbw; bounds.a.lower'];
    ubw             = [ubw; bounds.a.upper'];
    w0              = [w0;  guess.a(1,:)'];
    % Muscle-tendon forces
    FTtilde0        = MX.sym('FTtilde0',NMuscle);
    w               = [w {FTtilde0}];
    lbw             = [lbw; bounds.FTtilde.lower'];
    ubw             = [ubw; bounds.FTtilde.upper'];
    w0              = [w0;  guess.FTtilde(1,:)'];
    % Qs and Qdots
    X0              = MX.sym('X0',2*nq.all);
    w               = [w {X0}];    
    lbw             = [lbw; bounds.QsQdots.lower'];
    ubw             = [ubw; bounds.QsQdots.upper'];    
    w0              = [w0;  guess.QsQdots(1,:)'];
    % Arm activations
    a_a0            = MX.sym('a_a0',nq.arms);
    w               = [w {a_a0}];
    lbw             = [lbw; bounds.a_a.lower'];
    ubw             = [ubw; bounds.a_a.upper'];
    w0              = [w0;  guess.a_a(1,:)'];          
    % "Lift" initial conditions
    ak          = a0;
    FTtildek    = FTtilde0;
    Xk          = X0;
    a_ak        = a_a0;     
    % Time step
    h = (time_opt(2)-time_opt(1))/N;
    % Loop over mesh points
    for k=0:N-1
        % Define controls at mesh point (piecewise-constant in interval) 
        % Time derivative of muscle activations (states)
        vAk                 = MX.sym(['vA_' num2str(k)], NMuscle);
        w                   = [w {vAk}];
        lbw                 = [lbw; bounds.vA.lower'];
        ubw                 = [ubw; bounds.vA.upper'];
        w0                  = [w0; guess.vA(k+1,:)'];
        % Time derivative of muscle-tendon forces (states)
        dFTtildek           = MX.sym(['dFTtilde_' num2str(k)], NMuscle);
        w                   = [w {dFTtildek}];
        lbw                 = [lbw; bounds.dFTtilde.lower'];
        ubw                 = [ubw; bounds.dFTtilde.upper'];
        w0                  = [w0; guess.dFTtilde(k+1,:)'];  
        % Time derivative of Qdots (states) 
        Ak                  = MX.sym(['A_' num2str(k)], nq.all);
        w                   = [w {Ak}];
        lbw                 = [lbw; bounds.Qdotdots.lower'];
        ubw                 = [ubw; bounds.Qdotdots.upper'];
        w0                  = [w0; guess.Qdotdots(k+1,:)'];    
        % Arm excitations
        e_ak                = MX.sym(['e_a_' num2str(k)], nq.arms);
        w                   = [w {e_ak}];
        lbw                 = [lbw; bounds.e_a.lower'];
        ubw                 = [ubw; bounds.e_a.upper'];
        w0                  = [w0; guess.e_a(k+1,:)'];
        % Define states at collocation points    
        % Muscle activations
        akj = {};
        for j=1:d
            akj{j}  = MX.sym(['	a_' num2str(k) '_' num2str(j)], NMuscle);
            w       = {w{:}, akj{j}};
            lbw     = [lbw; bounds.a.lower'];
            ubw     = [ubw; bounds.a.upper'];
            w0      = [w0;  guess.a(k+1,:)'];
        end   
        % Muscle-tendon forces
        FTtildekj = {};
        for j=1:d
            FTtildekj{j} = ...
                MX.sym(['FTtilde_' num2str(k) '_' num2str(j)], NMuscle);
            w            = {w{:}, FTtildekj{j}};
            lbw          = [lbw; bounds.FTtilde.lower'];
            ubw          = [ubw; bounds.FTtilde.upper'];
            w0           = [w0;  guess.FTtilde(k+1,:)'];
        end
        % Qs and Qdots        
        Xkj = {};
        for j=1:d
            Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], 2*nq.all);
            w      = {w{:}, Xkj{j}};
            lbw    = [lbw; bounds.QsQdots.lower'];
            ubw    = [ubw; bounds.QsQdots.upper'];
            w0     = [w0;  guess.QsQdots(k+1,:)'];
        end   
        % Arm activations
        a_akj = {};
        for j=1:d
            a_akj{j}= MX.sym(['	a_a_' num2str(k) '_' num2str(j)], nq.arms);
            w       = {w{:}, a_akj{j}};
            lbw     = [lbw; bounds.a_a.lower'];
            ubw     = [ubw; bounds.a_a.upper'];
            w0      = [w0;  guess.a_a(k+1,:)'];
        end   
        % Unscale variables for later use
        Xk_nsc          = Xk.*scaling.QsQdots';
        FTtildek_nsc    = FTtildek.*(scaling.FTtilde');
        Ak_nsc          = Ak.*scaling.Qdotdots';
        for j=1:d
            Xkj_nsc{j} = Xkj{j}.*scaling.QsQdots';
            FTtildekj_nsc{j} = FTtildekj{j}.*scaling.FTtilde';
        end                  
        % Call external functions
        % The first external function (F1) returns linear and angular
        % velocities of the calcaneus and toes, positions of the origin of 
        % the calcaneus and toes, and transforms from the calcaneus and 
        % toes to the ground. These variables are used to compute the 
        % contact forces that will be input of the second external function 
        % (F2).
        [outk1] = F1(Xk_nsc); 
        % Organize the locations of the contact spheres
        locSphere_s1_l = [paramsCM_nsc(loci.s1.r.x),...
            loc.s1.y,-paramsCM_nsc(loci.s1.r.z)]';        
        locSphere_s1_r = [paramsCM_nsc(loci.s1.r.x),...
            loc.s1.y,paramsCM_nsc(loci.s1.r.z)]';        
        locSphere_s2_l = [paramsCM_nsc(loci.s2.r.x),...
            loc.s2.y,-paramsCM_nsc(loci.s2.r.z)]';
        locSphere_s2_r = [paramsCM_nsc(loci.s2.r.x),...
            loc.s2.y,paramsCM_nsc(loci.s2.r.z)]';        
        locSphere_s3_l = [paramsCM_nsc(loci.s3.r.x),...
            loc.s3.y,-paramsCM_nsc(loci.s3.r.z)]';
        locSphere_s3_r = [paramsCM_nsc(loci.s3.r.x),...
            loc.s3.y,paramsCM_nsc(loci.s3.r.z)]';        
        locSphere_s4_l = [paramsCM_nsc(loci.s4.r.x),...
            loc.s4.y,-paramsCM_nsc(loci.s4.r.z)]';
        locSphere_s4_r = [paramsCM_nsc(loci.s4.r.x),...
            loc.s4.y,paramsCM_nsc(loci.s4.r.z)]';          
        locSphere_s5_l = [paramsCM_nsc(loci.s5.r.x),...
            loc.s5.y,-paramsCM_nsc(loci.s5.r.z)]';
        locSphere_s5_r = [paramsCM_nsc(loci.s5.r.x),...
            loc.s5.y,paramsCM_nsc(loci.s5.r.z)]';          
        locSphere_s6_l = [paramsCM_nsc(loci.s6.r.x),...
            loc.s6.y,-paramsCM_nsc(loci.s6.r.z)]';
        locSphere_s6_r = [paramsCM_nsc(loci.s6.r.x),...
            loc.s6.y,paramsCM_nsc(loci.s6.r.z)]';          
        % Compute contact forces
        force_s1_l = f_contactForce(stif,paramsCM_nsc(radi.s1),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s1_l,...
            outk1(calcn.l.pos),outk1(calcn.l.v_lin),...
            outk1(calcn.l.omega),outk1(calcn.TR.R.l),outk1(calcn.TR.T.l));        
        force_s2_l = f_contactForce(stif,paramsCM_nsc(radi.s2),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s2_l,...
            outk1(calcn.l.pos),outk1(calcn.l.v_lin),...
            outk1(calcn.l.omega),outk1(calcn.TR.R.l),outk1(calcn.TR.T.l));        
        force_s3_l = f_contactForce(stif,paramsCM_nsc(radi.s3),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s3_l,...
            outk1(calcn.l.pos),outk1(calcn.l.v_lin),...
            outk1(calcn.l.omega),outk1(calcn.TR.R.l),outk1(calcn.TR.T.l));        
        force_s4_l = f_contactForce(stif,paramsCM_nsc(radi.s4),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s4_l,...
            outk1(toes.l.pos),outk1(toes.l.v_lin),...
            outk1(toes.l.omega),outk1(toes.TR.R.l),outk1(toes.TR.T.l));        
        force_s5_l = f_contactForce(stif,paramsCM_nsc(radi.s5),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s5_l,...
            outk1(calcn.l.pos),outk1(calcn.l.v_lin),...
            outk1(calcn.l.omega),outk1(calcn.TR.R.l),outk1(calcn.TR.T.l));        
        force_s6_l = f_contactForce(stif,paramsCM_nsc(radi.s6),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s6_l,...
            outk1(toes.l.pos),outk1(toes.l.v_lin),...
            outk1(toes.l.omega),outk1(toes.TR.R.l),outk1(toes.TR.T.l));        
        force_s1_r = f_contactForce(stif,paramsCM_nsc(radi.s1),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s1_r,...
            outk1(calcn.r.pos),outk1(calcn.r.v_lin),...
            outk1(calcn.r.omega),outk1(calcn.TR.R.r),outk1(calcn.TR.T.r));        
        force_s2_r = f_contactForce(stif,paramsCM_nsc(radi.s2),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s2_r,...
            outk1(calcn.r.pos),outk1(calcn.r.v_lin),...
            outk1(calcn.r.omega),outk1(calcn.TR.R.r),outk1(calcn.TR.T.r));        
        force_s3_r = f_contactForce(stif,paramsCM_nsc(radi.s3),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s3_r,...
            outk1(calcn.r.pos),outk1(calcn.r.v_lin),...
            outk1(calcn.r.omega),outk1(calcn.TR.R.r),outk1(calcn.TR.T.r));        
        force_s4_r = f_contactForce(stif,paramsCM_nsc(radi.s4),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s4_r,...
            outk1(toes.r.pos),outk1(toes.r.v_lin),...
            outk1(toes.r.omega),outk1(toes.TR.R.r),outk1(toes.TR.T.r));        
        force_s5_r = f_contactForce(stif,paramsCM_nsc(radi.s5),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s5_r,...
            outk1(calcn.r.pos),outk1(calcn.r.v_lin),...
            outk1(calcn.r.omega),outk1(calcn.TR.R.r),outk1(calcn.TR.T.r));        
        force_s6_r = f_contactForce(stif,paramsCM_nsc(radi.s6),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s6_r,...
            outk1(toes.r.pos),outk1(toes.r.v_lin),...
            outk1(toes.r.omega),outk1(toes.TR.R.r),outk1(toes.TR.T.r));
        in_F2 = [force_s1_l,force_s2_l,force_s3_l,force_s4_l,...
            force_s5_l,force_s6_l,force_s1_r,force_s2_r,force_s3_r,...
            force_s4_r,force_s5_r,force_s6_r,...       
            paramsCM_nsc(loci.s1.r.x),paramsCM_nsc(loci.s1.r.z),...
            paramsCM_nsc(loci.s2.r.x),paramsCM_nsc(loci.s2.r.z),...
            paramsCM_nsc(loci.s3.r.x),paramsCM_nsc(loci.s3.r.z),...
            paramsCM_nsc(loci.s4.r.x),paramsCM_nsc(loci.s4.r.z),...
            paramsCM_nsc(loci.s5.r.x),paramsCM_nsc(loci.s5.r.z),...
            paramsCM_nsc(loci.s6.r.x),paramsCM_nsc(loci.s6.r.z),...
            paramsCM_nsc(radi.s1:radi.s6)'];
        % The second external function (F2) returns joint torques, GRFs,
        % and GRMs based on joint states, controls, contact forces, and
        % several parameters of the contact models.
        [Tk] = F2([Xk_nsc;Ak_nsc;in_F2']);
        % Get muscle-tendon lengths, velocities, and moment arms
        % Left leg
        qin_l = [Xk_nsc(jointi.hip_flex.l*2-1,1),...
            Xk_nsc(jointi.hip_add.l*2-1,1), ...
            Xk_nsc(jointi.hip_rot.l*2-1,1), ...
            Xk_nsc(jointi.knee.l*2-1,1), ...
            Xk_nsc(jointi.ankle.l*2-1,1),...
            Xk_nsc(jointi.subt.l*2-1,1),...
            Xk_nsc(jointi.trunk.ext*2-1,1),...
            Xk_nsc(jointi.trunk.ben*2-1,1),...
            Xk_nsc(jointi.trunk.rot*2-1,1)];  
        qdotin_l = [Xk_nsc(jointi.hip_flex.l*2,1),...
            Xk_nsc(jointi.hip_add.l*2,1),...
            Xk_nsc(jointi.hip_rot.l*2,1),...
            Xk_nsc(jointi.knee.l*2,1),...
            Xk_nsc(jointi.ankle.l*2,1),...
            Xk_nsc(jointi.subt.l*2,1),...
            Xk_nsc(jointi.trunk.ext*2,1),...
            Xk_nsc(jointi.trunk.ben*2,1),...
            Xk_nsc(jointi.trunk.rot*2,1)];  
        [lMTk_l,vMTk_l,MA_l] = f_lMT_vMT_dM(qin_l,qdotin_l);    
        MA.hip_flex.l   =  MA_l(mai(1).mus.l',1);
        MA.hip_add.l    =  MA_l(mai(2).mus.l',2);
        MA.hip_rot.l    =  MA_l(mai(3).mus.l',3);
        MA.knee.l       =  MA_l(mai(4).mus.l',4);
        MA.ankle.l      =  MA_l(mai(5).mus.l',5);  
        MA.subt.l       =  MA_l(mai(6).mus.l',6); 
        % For the back muscles, we want left and right together: left
        % first, right second. In MuscleInfo, we first have the right
        % muscles (44:46) and then the left muscles (47:49). Since the back
        % muscles only depend on back dofs, we do not care if we extract
        % them "from the left or right leg" so here we just picked left.
        MA.trunk_ext    =  MA_l([47:49,mai(7).mus.l]',7);
        MA.trunk_ben    =  MA_l([47:49,mai(8).mus.l]',8);
        MA.trunk_rot    =  MA_l([47:49,mai(9).mus.l]',9);
        % Right leg
        qin_r = [Xk_nsc(jointi.hip_flex.r*2-1,1),...
            Xk_nsc(jointi.hip_add.r*2-1,1),...
            Xk_nsc(jointi.hip_rot.r*2-1,1),...
            Xk_nsc(jointi.knee.r*2-1,1),...
            Xk_nsc(jointi.ankle.r*2-1,1),...
            Xk_nsc(jointi.subt.r*2-1,1),...
            Xk_nsc(jointi.trunk.ext*2-1,1),...
            Xk_nsc(jointi.trunk.ben*2-1,1),...
            Xk_nsc(jointi.trunk.rot*2-1,1)];  
        qdotin_r = [Xk_nsc(jointi.hip_flex.r*2,1),...
            Xk_nsc(jointi.hip_add.r*2,1),...
            Xk_nsc(jointi.hip_rot.r*2,1),...
            Xk_nsc(jointi.knee.r*2,1),...
            Xk_nsc(jointi.ankle.r*2,1),...
            Xk_nsc(jointi.subt.r*2,1),...
            Xk_nsc(jointi.trunk.ext*2,1),...
            Xk_nsc(jointi.trunk.ben*2,1),...
            Xk_nsc(jointi.trunk.rot*2,1)];      
        [lMTk_r,vMTk_r,MA_r] = f_lMT_vMT_dM(qin_r,qdotin_r);
        % Here we take the indices from left since the vector is 1:49
        MA.hip_flex.r   =  MA_r(mai(1).mus.l',1);
        MA.hip_add.r    =  MA_r(mai(2).mus.l',2);
        MA.hip_rot.r    =  MA_r(mai(3).mus.l',3);
        MA.knee.r       =  MA_r(mai(4).mus.l',4);
        MA.ankle.r      =  MA_r(mai(5).mus.l',5);
        MA.subt.r       =  MA_r(mai(6).mus.l',6);
        % Both legs
        % In MuscleInfo, we first have the right back muscles (44:46) and 
        % then the left back muscles (47:49). Here we re-organize so that
        % we have first the left muscles and then the right muscles.
        lMTk_lr     = [lMTk_l([1:43,47:49],1);lMTk_r(1:46,1)];
        vMTk_lr     = [vMTk_l([1:43,47:49],1);vMTk_r(1:46,1)];   
        % Get muscle-tendon forces and derive Hill-equilibrium       
        [Hilldiffk,FTk,~,~,~,~] =  f_forceEquilibrium_FtildeState_all(...
                ak,FTtildek.*scaling.FTtilde',...
                dFTtildek.*scaling.dFTtilde,lMTk_lr,vMTk_lr,tensions);
        % Get passive torques
        Tau_passk.hip.flex.l    = f_PassiveMoments(k_pass.hip.flex,...
            theta.pass.hip.flex,Xk_nsc(jointi.hip_flex.l*2-1,1),...
            Xk_nsc(jointi.hip_flex.l*2,1));
        Tau_passk.hip.flex.r    = f_PassiveMoments(k_pass.hip.flex,...
            theta.pass.hip.flex,Xk_nsc(jointi.hip_flex.r*2-1,1),...
            Xk_nsc(jointi.hip_flex.r*2,1));
        Tau_passk.hip.add.l     = f_PassiveMoments(k_pass.hip.add,...
            theta.pass.hip.add,Xk_nsc(jointi.hip_add.l*2-1,1),...
            Xk_nsc(jointi.hip_add.l*2,1));
        Tau_passk.hip.add.r     = f_PassiveMoments(k_pass.hip.add,...
            theta.pass.hip.add,Xk_nsc(jointi.hip_add.r*2-1,1),...
            Xk_nsc(jointi.hip_add.r*2,1));
        Tau_passk.hip.rot.l     = f_PassiveMoments(k_pass.hip.rot,...
            theta.pass.hip.rot,Xk_nsc(jointi.hip_rot.l*2-1,1),...
            Xk_nsc(jointi.hip_rot.l*2,1));
        Tau_passk.hip.rot.r     = f_PassiveMoments(k_pass.hip.rot,...
            theta.pass.hip.rot,Xk_nsc(jointi.hip_rot.r*2-1,1),...
            Xk_nsc(jointi.hip_rot.r*2,1));
        Tau_passk.knee.l        = f_PassiveMoments(k_pass.knee,...
            theta.pass.knee,Xk_nsc(jointi.knee.l*2-1,1),...
            Xk_nsc(jointi.knee.l*2,1));
        Tau_passk.knee.r        = f_PassiveMoments(k_pass.knee,...
            theta.pass.knee,Xk_nsc(jointi.knee.r*2-1,1),...
            Xk_nsc(jointi.knee.r*2,1));
        Tau_passk.ankle.l       = f_PassiveMoments(k_pass.ankle,...
            theta.pass.ankle,Xk_nsc(jointi.ankle.l*2-1,1),...
            Xk_nsc(jointi.ankle.l*2,1));
        Tau_passk.ankle.r       = f_PassiveMoments(k_pass.ankle,...
            theta.pass.ankle,Xk_nsc(jointi.ankle.r*2-1,1),...
            Xk_nsc(jointi.ankle.r*2,1));        
        Tau_passk.subt.l       = f_PassiveMoments(k_pass.subt,...
            theta.pass.subt,Xk_nsc(jointi.subt.l*2-1,1),...
            Xk_nsc(jointi.subt.l*2,1));
        Tau_passk.subt.r       = f_PassiveMoments(k_pass.subt,...
            theta.pass.subt,Xk_nsc(jointi.subt.r*2-1,1),...
            Xk_nsc(jointi.subt.r*2,1));        
        Tau_passk.trunk.ext     = f_PassiveMoments(k_pass.trunk.ext,...
            theta.pass.trunk.ext,Xk_nsc(jointi.trunk.ext*2-1,1),...
            Xk_nsc(jointi.trunk.ext*2,1));
        Tau_passk.trunk.ben     = f_PassiveMoments(k_pass.trunk.ben,...
            theta.pass.trunk.ben,Xk_nsc(jointi.trunk.ben*2-1,1),...
            Xk_nsc(jointi.trunk.ben*2,1));
        Tau_passk.trunk.rot     = f_PassiveMoments(k_pass.trunk.rot,...
            theta.pass.trunk.rot,Xk_nsc(jointi.trunk.rot*2-1,1),...
            Xk_nsc(jointi.trunk.rot*2,1));        
        Tau_passk_all = [Tau_passk.hip.flex.l,Tau_passk.hip.flex.r,...
            Tau_passk.hip.add.l,Tau_passk.hip.add.r,...
            Tau_passk.hip.rot.l,Tau_passk.hip.rot.r,...
            Tau_passk.knee.l,Tau_passk.knee.r,Tau_passk.ankle.l,...
            Tau_passk.ankle.r,Tau_passk.subt.l,Tau_passk.subt.r,...
            Tau_passk.trunk.ext,Tau_passk.trunk.ben,...
            Tau_passk.trunk.rot]';
        % Loop over collocation points
        Xk_nsc_end          = D(1)*Xk_nsc;
        FTtildek_nsc_end    = D(1)*FTtildek_nsc;
        ak_end              = D(1)*ak;
        a_ak_end            = D(1)*a_ak;
        for j=1:d
            % Expression for the state derivatives at the collocation point
            xp_nsc          = C(1,j+1)*Xk_nsc;
            FTtildep_nsc    = C(1,j+1)*FTtildek_nsc;
            ap              = C(1,j+1)*ak;
            a_ap            = C(1,j+1)*a_ak;
            for r=1:d
                xp_nsc       = xp_nsc + C(r+1,j+1)*Xkj_nsc{r};
                FTtildep_nsc = FTtildep_nsc + C(r+1,j+1)*FTtildekj_nsc{r};
                ap           = ap + C(r+1,j+1)*akj{r};
                a_ap         = a_ap + C(r+1,j+1)*a_akj{r};
            end 
            % Append collocation equations
            % Dynamic constraints are scaled using the same scale
            % factors as was used to scale the states
            % Activation dynamics (implicit formulation)  
            g       = {g{:}, (h*vAk.*scaling.vA - ap)./scaling.a};
            lbg     = [lbg; zeros(NMuscle,1)];
            ubg     = [ubg; zeros(NMuscle,1)]; 
            % Contraction dynamics (implicit formulation)               
            g       = {g{:}, (h*dFTtildek.*scaling.dFTtilde - ...
                FTtildep_nsc)./(scaling.FTtilde')};
            lbg     = [lbg; zeros(NMuscle,1)];
            ubg     = [ubg; zeros(NMuscle,1)];
            % Skeleton dynamics (implicit formulation)        
            xj_nsc  = [...
                Xkj_nsc{j}(2); Ak_nsc(1); Xkj_nsc{j}(4); Ak_nsc(2);...
                Xkj_nsc{j}(6); Ak_nsc(3); Xkj_nsc{j}(8); Ak_nsc(4);...
                Xkj_nsc{j}(10); Ak_nsc(5); Xkj_nsc{j}(12); Ak_nsc(6);...
                Xkj_nsc{j}(14); Ak_nsc(7); Xkj_nsc{j}(16); Ak_nsc(8);...
                Xkj_nsc{j}(18); Ak_nsc(9); Xkj_nsc{j}(20); Ak_nsc(10);...
                Xkj_nsc{j}(22); Ak_nsc(11); Xkj_nsc{j}(24); Ak_nsc(12);...
                Xkj_nsc{j}(26); Ak_nsc(13); Xkj_nsc{j}(28); Ak_nsc(14);...
                Xkj_nsc{j}(30); Ak_nsc(15); Xkj_nsc{j}(32); Ak_nsc(16);...
                Xkj_nsc{j}(34); Ak_nsc(17); Xkj_nsc{j}(36); Ak_nsc(18);...
                Xkj_nsc{j}(38); Ak_nsc(19); Xkj_nsc{j}(40); Ak_nsc(20);...
                Xkj_nsc{j}(42); Ak_nsc(21); Xkj_nsc{j}(44); Ak_nsc(22);...
                Xkj_nsc{j}(46); Ak_nsc(23); Xkj_nsc{j}(48); Ak_nsc(24); ...
                Xkj_nsc{j}(50); Ak_nsc(25); Xkj_nsc{j}(52); Ak_nsc(26); ...
                Xkj_nsc{j}(54); Ak_nsc(27); Xkj_nsc{j}(56); Ak_nsc(28); ...
                Xkj_nsc{j}(58); Ak_nsc(29);];
            g       = {g{:}, (h*xj_nsc - xp_nsc)./(scaling.QsQdots')};
            lbg     = [lbg; zeros(2*nq.all,1)];
            ubg     = [ubg; zeros(2*nq.all,1)];   
            % Arm activation dynamics (explicit formulation)   
            dadt    = f_ArmActivationDynamics(e_ak,a_akj{j});
            g       = {g{:}, (h*dadt - a_ap)./scaling.a_a};
            lbg     = [lbg; zeros(nq.arms,1)];
            ubg     = [ubg; zeros(nq.arms,1)]; 
            % Add contribution to the end state
            Xk_nsc_end = Xk_nsc_end + D(j+1)*Xkj_nsc{j};
            FTtildek_nsc_end = FTtildek_nsc_end + D(j+1)*FTtildekj_nsc{j};
            ak_end = ak_end + D(j+1)*akj{j};  
            a_ak_end = a_ak_end + D(j+1)*a_akj{j};    
            % Add contribution to quadrature function
            J = J + ...
                W.Qs*B(j+1)*(f_J28(Xk(Qsi(residual_bptyi),1)-...
                    Qs.allinterpfilt(k+1,residual_bptyi+1)'... 
                    ./scaling.Qs(residual_bptyi)'))*h +...
                W.GRF*B(j+1)*(f_J6((Tk(GRFi.all,1)./scaling.GRF')-...
                    GRF.val.allinterp(k+1,2:end)'./scaling.GRF'))*h +...
                W.GRM*B(j+1)*(f_J6((Tk(GRMi.all,1)./scaling.GRM')-...
                    GRF.MorGF.allinterp(k+1,2:end)'./scaling.GRM'))*h +...
                W.ID_act*B(j+1)*(f_J23((Tk(residuals_acti,1)./scaling.T(1)')-...
                    ID.allinterp(k+1,2+nq.abs:end)'./scaling.T(1)))*h +...
                W.a*B(j+1)*(f_J92(akj{j}))*h + ...
                W.u*B(j+1)*(f_J29(Ak))*h +...
                W.u*B(j+1)*(f_J92(vAk))*h +...
                W.u*B(j+1)*(f_J92(dFTtildek))*h;
        end                              
        % Add path constraints
        % Pelvis residuals (same as from inverse dynamics)
        g   = {g{:},(ID.allinterp(k+1,2:7)' - ...
            Tk(ground_pelvisi,1))./scaling.T(1)};
        lbg = [lbg; zeros(nq.abs,1)];
        ubg = [ubg; zeros(nq.abs,1)];    
        % Muscle-driven joint torques for the lower limbs and the trunk
        % Hip flexion, left
        Ft_hip_flex_l   = FTk(mai(1).mus.l',1);
        T_hip_flex_l    = f_T27(MA.hip_flex.l,Ft_hip_flex_l);
        g               = {g{:},Tk(jointi.hip_flex.l,1) - ...
            (T_hip_flex_l + Tau_passk.hip.flex.l)};        
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];    
        % Hip flexion, right
        Ft_hip_flex_r   = FTk(mai(1).mus.r',1);
        T_hip_flex_r    = f_T27(MA.hip_flex.r,Ft_hip_flex_r);
        g               = {g{:},Tk(jointi.hip_flex.r,1) - ...
            (T_hip_flex_r + Tau_passk.hip.flex.r)};        
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];   
        % Hip adduction, left
        Ft_hip_add_l    = FTk(mai(2).mus.l',1);
        T_hip_add_l     = f_T27(MA.hip_add.l,Ft_hip_add_l);
        g               = {g{:},Tk(jointi.hip_add.l,1) - ...
            (T_hip_add_l + Tau_passk.hip.add.l)};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];    
        % Hip adduction, right
        Ft_hip_add_r    = FTk(mai(2).mus.r',1);
        T_hip_add_r     = f_T27(MA.hip_add.r,Ft_hip_add_r);
        g               = {g{:},Tk(jointi.hip_add.r,1) - ...
            (T_hip_add_r + Tau_passk.hip.add.r)};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];  
        % Hip rotation, left
        Ft_hip_rot_l    = FTk(mai(3).mus.l',1);
        T_hip_rot_l     = f_T27(MA.hip_rot.l,Ft_hip_rot_l);
        g               = {g{:},Tk(jointi.hip_rot.l,1) - ...
            (T_hip_rot_l + Tau_passk.hip.rot.l)};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];    
        % Hip rotation, right
        Ft_hip_rot_r    = FTk(mai(3).mus.r',1);
        T_hip_rot_r     = f_T27(MA.hip_rot.r,Ft_hip_rot_r);
        g               = {g{:},Tk(jointi.hip_rot.r,1) - ...
            (T_hip_rot_r + Tau_passk.hip.rot.r)};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];   
        % Knee, left
        Ft_knee_l       = FTk(mai(4).mus.l',1);
        T_knee_l        = f_T13(MA.knee.l,Ft_knee_l);
        g               = {g{:},Tk(jointi.knee.l,1) - ...
            (T_knee_l + Tau_passk.knee.l)};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];    
        % Knee, right
        Ft_knee_r       = FTk(mai(4).mus.r',1);
        T_knee_r        = f_T13(MA.knee.r,Ft_knee_r);
        g               = {g{:},Tk(jointi.knee.r,1) - ...
            (T_knee_r + Tau_passk.knee.r)};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];    
        % Ankle, left
        Ft_ankle_l      = FTk(mai(5).mus.l',1);
        T_ankle_l       = f_T12(MA.ankle.l,Ft_ankle_l);
        g               = {g{:},Tk(jointi.ankle.l,1) - ...
            (T_ankle_l + Tau_passk.ankle.l)};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];    
        % Ankle, right
        Ft_ankle_r      = FTk(mai(5).mus.r',1);
        T_ankle_r       = f_T12(MA.ankle.r,Ft_ankle_r);
        g               = {g{:},Tk(jointi.ankle.r,1) - ...
            (T_ankle_r + Tau_passk.ankle.r)};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];
        % Subtalar, left
        Ft_subt_l       = FTk(mai(6).mus.l',1);
        T_subt_l        = f_T12(MA.subt.l,Ft_subt_l);
        g               = {g{:},(Tk(jointi.subt.l,1) - ...
            (T_subt_l + Tau_passk.subt.l))};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];    
        % Subtalar, right
        Ft_subt_r       = FTk(mai(6).mus.r',1);
        T_subt_r        = f_T12(MA.subt.r,Ft_subt_r);
        g               = {g{:},(Tk(jointi.subt.r,1) - ...
            (T_subt_r + Tau_passk.subt.r))};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];
        % Lumbar extension
        Ft_trunk_ext    = FTk([mai(7).mus.l,mai(7).mus.r]',1);
        T_trunk_ext     = f_T6(MA.trunk_ext,Ft_trunk_ext);
        g               = {g{:},Tk(jointi.trunk.ext,1) - ...
            (T_trunk_ext + Tau_passk.trunk.ext)};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];
        % Lumbar bending
        Ft_trunk_ben    = FTk([mai(8).mus.l,mai(8).mus.r]',1);
        T_trunk_ben     = f_T6(MA.trunk_ben,Ft_trunk_ben);
        g               = {g{:},Tk(jointi.trunk.ben,1) - ...
            (T_trunk_ben + Tau_passk.trunk.ben)};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];
        % Lumbar rotation
        Ft_trunk_rot    = FTk([mai(9).mus.l,mai(9).mus.r]',1);
        T_trunk_rot     = f_T6(MA.trunk_rot,Ft_trunk_rot);
        g               = {g{:},Tk(jointi.trunk.rot,1) - ...
            (T_trunk_rot + Tau_passk.trunk.rot)};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];
        % Torque-driven joint torques for the arms
        % Arms
        g               = {g{:},Tk(armsi,1)/scaling.ArmTau - a_ak};
        lbg             = [lbg; zeros(nq.arms,1)];
        ubg             = [ubg; zeros(nq.arms,1)];
        % Activation dynamics (implicit formulation)
        tact = 0.015;
        tdeact = 0.06;
        act1 = vAk*scaling.vA + ak./(ones(size(ak,1),1)*tdeact);
        act2 = vAk*scaling.vA + ak./(ones(size(ak,1),1)*tact);
        % act1
        g               = {g{:},act1};
        lbg             = [lbg; zeros(NMuscle,1)];
        ubg             = [ubg; inf*ones(NMuscle,1)]; 
        % act2
        g               = {g{:},act2};
        lbg             = [lbg; -inf*ones(NMuscle,1)];
        ubg             = [ubg; ones(NMuscle,1)./(ones(NMuscle,1)*tact)];        
        % Contraction dynamics (implicit formulation)
        g               = {g{:},Hilldiffk};
        lbg             = [lbg; zeros(NMuscle,1)];
        ubg             = [ubg; zeros(NMuscle,1)];      
        % New NLP variables for states at end of interval
        if k ~= N-1
            % Muscle activations
            ak              = MX.sym(['a_' num2str(k+1)], NMuscle);
            w               = {w{:}, ak};
            lbw             = [lbw; bounds.a.lower'];
            ubw             = [ubw; bounds.a.upper'];
            w0              = [w0;  guess.a(k+2,:)'];
            % Muscle-tendon forces
            FTtildek        = MX.sym(['FTtilde_' num2str(k+1)], NMuscle);
            w               = {w{:}, FTtildek};
            lbw             = [lbw; bounds.FTtilde.lower'];
            ubw             = [ubw; bounds.FTtilde.upper'];
            w0              = [w0;  guess.FTtilde(k+2,:)'];    
            % Qs and Qdots
            Xk              = MX.sym(['X_' num2str(k+1)],2*nq.all);
            w               = {w{:}, Xk};
            lbw             = [lbw; bounds.QsQdots.lower'];
            ubw             = [ubw; bounds.QsQdots.upper']; 
            w0              = [w0;  guess.QsQdots(k+2,:)'];
            % Arm activations
            a_ak            = MX.sym(['a_a_' num2str(k+1)], nq.arms);
            w               = {w{:}, a_ak};
            lbw             = [lbw; bounds.a_a.lower'];
            ubw             = [ubw; bounds.a_a.upper'];
            w0              = [w0;  guess.a_a(k+2,:)'];
        else
            % Muscle activations
            ak              = MX.sym(['a_' num2str(k+1)], NMuscle);
            w               = {w{:}, ak};
            lbw             = [lbw; bounds.a.lower'];
            ubw             = [ubw; bounds.a.upper'];
            w0              = [w0;  guess.a(end,:)'];
            % Muscle-tendon forces
            FTtildek        = MX.sym(['FTtilde_' num2str(k+1)], NMuscle);
            w               = {w{:}, FTtildek};
            lbw             = [lbw; bounds.FTtilde.lower'];
            ubw             = [ubw; bounds.FTtilde.upper'];
            w0              = [w0;  guess.FTtilde(end,:)'];    
            % Qs and Qdots
            Xk              = MX.sym(['X_' num2str(k+1)],2*nq.all);
            w               = {w{:}, Xk};
            lbw             = [lbw; bounds.QsQdots.lower'];
            ubw             = [ubw; bounds.QsQdots.upper']; 
            w0              = [w0;  guess.QsQdots(end,:)'];
            % Arm activations
            a_ak            = MX.sym(['a_a_' num2str(k+1)], nq.arms);
            w               = {w{:}, a_ak};
            lbw             = [lbw; bounds.a_a.lower'];
            ubw             = [ubw; bounds.a_a.upper'];
            w0              = [w0;  guess.a_a(end,:)'];
        end
        % Rescale variables to impose equality constraints
        Xk_end = (Xk_nsc_end)./scaling.QsQdots';
        FTtildek_end = (FTtildek_nsc_end)./scaling.FTtilde';
        % Add equality constraints (next interval starts with end values of 
        % states from previous interval)
        g   = {g{:}, Xk_end-Xk, FTtildek_end-FTtildek, ...
            ak_end-ak, a_ak_end-a_ak};
        lbg = [lbg; zeros(2*nq.all + NMuscle + NMuscle + nq.arms,1)];
        ubg = [ubg; zeros(2*nq.all + NMuscle + NMuscle + nq.arms,1)];    

    end   
    
    % Assert bounds / IG
    % Lower bounds smaller than upper bounds
    assert_bw = isempty(find(lbw <= ubw == 0,1));
    assert_bg = isempty(find(lbg <= ubg == 0,1));
    % Design variables between -1 and 1    
    assert_bwl = isempty(find(lbw < -1 == 1,1));
    assert_bwu = isempty(find(1 < ubw == 1,1));
    % Initial guess within bounds
    assert_w0_ubw = isempty(find(w0 <= ubw == 0,1));
    assert_w0_lbw = isempty(find(lbw <= w0 == 0,1)); 
                
    % Create an NLP solver
    prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy      = 'adaptive';
    options.ipopt.max_iter = 10000;
    options.ipopt.tol = 1*10^(-tol_ipopt);
    solver = nlpsol('solver', 'ipopt', prob, options);
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
    % Solve problem
    sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
        'lbg', lbg, 'ubg', ubg);    
    diary off
    % Extract results
    w_opt = full(sol.x);
    g_opt = full(sol.g);  
    % Save results
    setup.tolerance.ipopt = tol_ipopt;
    setup.bounds = bounds;
    setup.scaling = scaling;
    setup.guess = guess;
    setup.lbw = lbw;
    setup.ubw = ubw;
    % Save results and setup
    save([pathresults,'/',namescript,'/w',savename],'w_opt');
    save([pathresults,'/',namescript,'/g',savename],'g_opt');
    save([pathresults,'/',namescript,'/s',savename],'setup');
end

%% Analyze results
if analyseResults
    %% Load results
    if loadResults
        p = mfilename('fullpath');
        [~,namescript,~] = fileparts(p);
        pathresults = [pathRepo,'/Results'];
        load([pathresults,'/',namescript,'/w',savename]);
        load([pathresults,'/',namescript,'/g',savename]);
        load([pathresults,'/',namescript,'/s',savename]);
    end
    
    %% Extract results
    % All optimized design variables are saved in a single column vector      
    % Number of design variables    
    NControls = NMuscle+NMuscle+nq.all+nq.arms;
    NStates = NMuscle+NMuscle+2*nq.all+nq.arms;
    NParameters = np;
    % In the loop
    Nwl = NControls+d*NStates+NStates;
    % In total
    Nw = NParameters+NStates+N*Nwl;
    % Before the variable corresponding to the first collocation point
    Nwm = NParameters+NStates+NControls;
    % Here we extract the results and re-organize them for analysis  
    % Static parameters
    paramsCM_opt    = w_opt(1:NParameters);  
    % Mesh points
    % Muscle activations and muscle-tendon forces
    a_opt = zeros(N+1,NMuscle);
    FTtilde_opt = zeros(N+1,NMuscle);
    for i = 1:NMuscle
        a_opt(:,i) = w_opt(NParameters+i:Nwl:Nw);
        FTtilde_opt(:,i) = w_opt(NParameters+NMuscle+i:Nwl:Nw);
    end
    % Qs and Qdots
    q_opt = zeros(N+1,nq.all);
    qdot_opt = zeros(N+1,nq.all);
    count = 0;
    for i = 1:2:2*nq.all
        count = count +1;
        q_opt(:,count)      = w_opt(NParameters+NMuscle+NMuscle+i:Nwl:Nw);
        qdot_opt(:,count)   = w_opt(NParameters+NMuscle+NMuscle+i+1:Nwl:Nw);
    end
    % Arm activations
    a_a_opt = zeros(N+1,nq.arms);
    for i = 1:nq.arms
        a_a_opt(:,i) = w_opt(NParameters+NMuscle+NMuscle+2*nq.all+i:Nwl:Nw);
    end
    % Time derivative of muscle activations and muscle-tendon forces
    vA_opt = zeros(N,NMuscle);
    dFTtilde_opt = zeros(N,NMuscle);
    for i = 1:NMuscle
        vA_opt(:,i)         = w_opt(NParameters+NMuscle+NMuscle+2*nq.all+...
            nq.arms+i:Nwl:Nw);
        dFTtilde_opt(:,i)   = w_opt(NParameters+NMuscle+NMuscle+2*nq.all+...
            nq.arms+NMuscle+i:Nwl:Nw);
    end
    % Time derivative of joint velocities
    qdotdot_opt = zeros(N,nq.all);
    for i = 1:nq.all
        qdotdot_opt(:,i)    = w_opt(NParameters+NMuscle+NMuscle+2*nq.all+...
            nq.arms+NMuscle+NMuscle+i:Nwl:Nw);
    end
    % Arm excitations
    e_a_opt = zeros(N,nq.arms);
    for i = 1:nq.arms
        e_a_opt(:,i)       = w_opt(NParameters+NMuscle+NMuscle+2*nq.all+...
            nq.arms+NMuscle+NMuscle+nq.all+i:Nwl:Nw); 
    end    
    % Collocation points
    % Muscle activations
    a_opt_ext=zeros(N*(d+1)+1,NMuscle);
    a_opt_ext(1:(d+1):end,:)= a_opt;
    for nmusi=1:NMuscle
        a_opt_ext(2:(d+1):end,nmusi) = w_opt(Nwm+nmusi:Nwl:Nw);
        a_opt_ext(3:(d+1):end,nmusi) = ...
            w_opt(Nwm+NMuscle+nmusi:Nwl:Nw);
        a_opt_ext(4:(d+1):end,nmusi) = ...
            w_opt(Nwm+NMuscle+NMuscle+nmusi:Nwl:Nw);
    end
    % Muscle-tendon forces
    FTtilde_opt_ext=zeros(N*(d+1)+1,NMuscle);
    FTtilde_opt_ext(1:(d+1):end,:)= FTtilde_opt;
    for nmusi=1:NMuscle
        FTtilde_opt_ext(2:(d+1):end,nmusi) = ...
            w_opt(Nwm+d*NMuscle+nmusi:Nwl:Nw);
        FTtilde_opt_ext(3:(d+1):end,nmusi) = ...
            w_opt(Nwm+d*NMuscle+NMuscle+nmusi:Nwl:Nw);
        FTtilde_opt_ext(4:(d+1):end,nmusi) = ...
            w_opt(Nwm+d*NMuscle+NMuscle+NMuscle+nmusi:Nwl:Nw);
    end
    % Qs and Qdots
    q_opt_ext=zeros(N*(d+1)+1,nq.all);
    q_opt_ext(1:(d+1):end,:)= q_opt;
    q_dot_opt_ext=zeros(N*(d+1)+1,nq.all);
    q_dot_opt_ext(1:(d+1):end,:)= qdot_opt;
    nqi_col = 1:2:2*nq.all;
    for nqi=1:nq.all
        nqi_q = nqi_col(nqi);
        q_opt_ext(2:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+nqi_q:Nwl:Nw);   
        q_opt_ext(3:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+2*nq.all+nqi_q:Nwl:Nw);  
        q_opt_ext(4:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+2*nq.all+2*nq.all+nqi_q:Nwl:Nw);  
        q_dot_opt_ext(2:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+nqi_q+1:Nwl:Nw);   
        q_dot_opt_ext(3:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+2*nq.all+nqi_q+1:Nwl:Nw);  
        q_dot_opt_ext(4:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+2*nq.all+2*nq.all+nqi_q+1:Nwl:Nw);
    end
    % Arm activations
    a_a_opt_ext=zeros(N*(d+1)+1,nq.arms);
    a_a_opt_ext(1:(d+1):end,:)= a_a_opt;
    for nmusi=1:nq.arms
        a_a_opt_ext(2:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+d*2*nq.all+nmusi:Nwl:Nw);
        a_a_opt_ext(3:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+d*2*nq.all+nq.arms+nmusi:Nwl:Nw);
        a_a_opt_ext(4:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+d*2*nq.all+nq.arms+nq.arms+nmusi:Nwl:Nw);
    end
    
    %% Unscale results
    % Parameters
    paramsCM_opt_nsc = f_nsc18(paramsCM_opt,scaling.params.v,scaling.params.r);
    % States at mesh points
    % Qs (1:N-1)
    q_opt_unsc.rad = q_opt(1:end-1,:).*repmat(...
    scaling.Qs,size(q_opt(1:end-1,:),1),1); 
    % Convert in degrees
    q_opt_unsc.deg = q_opt_unsc.rad;
    q_opt_unsc.deg(:,[1:3,7:end]) = q_opt_unsc.deg(:,[1:3,7:end]).*180/pi;
    % Qs (1:N)
    q_opt_unsc_all.rad = q_opt(1:end,:).*repmat(...
        scaling.Qs,size(q_opt(1:end,:),1),1); 
    % Convert in degrees
    q_opt_unsc_all.deg = q_opt_unsc_all.rad;
    q_opt_unsc_all.deg(:,[1:3,7:end]) = ...
        q_opt_unsc_all.deg(:,[1:3,7:end]).*180/pi;
    % Qdots (1:N-1)
    qdot_opt_unsc.rad = qdot_opt(1:end-1,:).*repmat(...
        scaling.Qdots,size(qdot_opt(1:end-1,:),1),1);
    % Convert in degrees
    qdot_opt_unsc.deg = qdot_opt_unsc.rad;
    qdot_opt_unsc.deg(:,[1:3,7:end]) = ...
        qdot_opt_unsc.deg(:,[1:3,7:end]).*180/pi;
    % Qdots (1:N)
    qdot_opt_unsc_all.rad = qdot_opt(1:end,:).*repmat(...
        scaling.Qdots,size(qdot_opt(1:end,:),1),1); 
    % Muscle activations
    a_opt_unsc = a_opt(1:end-1,:).*repmat(...
        scaling.a,size(a_opt(1:end-1,:),1),size(a_opt,2));
    % Muscle-tendon forces
    FTtilde_opt_unsc = FTtilde_opt(1:end-1,:).*repmat(...
        scaling.FTtilde,size(FTtilde_opt(1:end-1,:),1),1);
    % Arm activations
    a_a_opt_unsc = a_a_opt(1:end-1,:).*repmat(...
        scaling.a_a,size(a_a_opt(1:end-1,:),1),size(a_a_opt,2));
    % Controls at mesh points
    % Time derivative of Qdots
    qdotdot_opt_unsc.rad = ...
        qdotdot_opt.*repmat(scaling.Qdotdots,size(qdotdot_opt,1),1);
    % Convert in degrees
    qdotdot_opt_unsc.deg = qdotdot_opt_unsc.rad;
    qdotdot_opt_unsc.deg(:,[1:3,7:end]) = ...
        qdotdot_opt_unsc.deg(:,[1:3,7:end]).*180/pi;
    % Time derivative of muscle activations (states)
    vA_opt_unsc = vA_opt.*repmat(scaling.vA,size(vA_opt,1),size(vA_opt,2));
    tact = 0.015;
    tdeact = 0.06;
    % Get muscle excitations from time derivative of muscle activations
    e_opt_unsc = computeExcitationRaasch(a_opt_unsc,vA_opt_unsc,...
        ones(1,NMuscle)*tdeact,ones(1,NMuscle)*tact);
    % Time derivative of muscle-tendon forces
    dFTtilde_opt_unsc = dFTtilde_opt.*repmat(...
        scaling.dFTtilde,size(dFTtilde_opt,1),size(dFTtilde_opt,2));
    % Arm excitations
    e_a_opt_unsc = e_a_opt.*repmat(scaling.e_a,size(e_a_opt,1),...
        size(e_a_opt,2));
    
    %% Time grid    
    % Mesh points
    tgrid = linspace(time_opt(1),time_opt(end),N+1);
    dtime = zeros(1,d+1);
    for i=1:4
        dtime(i)=tau_root(i)*((time_opt(end)-time_opt(1))/N);
    end
    % Mesh points and collocation points
    tgrid_ext = zeros(1,(d+1)*N+1);
    for i=1:N
        tgrid_ext(((i-1)*4+1):1:i*4)=tgrid(i)+dtime;
    end
    tgrid_ext(end)=time_opt(end); 
 
    %% Joint torques, ground reaction forces and moments at opt solution
    Xk_Qs_Qdots_opt             = zeros(N,2*size(q_opt_unsc.rad,2));
    Xk_Qs_Qdots_opt(:,1:2:end)  = q_opt_unsc.rad;
    Xk_Qs_Qdots_opt(:,2:2:end)  = qdot_opt_unsc.rad;
    Xk_Qdotdots_opt             = qdotdot_opt_unsc.rad;  
    out1_res_opt                = zeros(N,toes.TR.T.r(end));
    out2_res_opt                = zeros(N,GRMi.l(end));
    for i = 1:N
        out1_res = F1(Xk_Qs_Qdots_opt(i,:));
        out1_res_opt(i,:) = full(out1_res);
        % Re-organise the contact sphere locations
        locSphere_s1_l_opt = [paramsCM_opt_nsc(loci.s1.r.x),...
            loc.s1.y,-paramsCM_opt_nsc(loci.s1.r.z)]';
        locSphere_s1_r_opt = [paramsCM_opt_nsc(loci.s1.r.x),...
            loc.s1.y,paramsCM_opt_nsc(loci.s1.r.z)]';
        locSphere_s2_l_opt = [paramsCM_opt_nsc(loci.s2.r.x),...
            loc.s2.y,-paramsCM_opt_nsc(loci.s2.r.z)]';
        locSphere_s2_r_opt = [paramsCM_opt_nsc(loci.s2.r.x),...
            loc.s2.y,paramsCM_opt_nsc(loci.s2.r.z)]';
        locSphere_s3_l_opt = [paramsCM_opt_nsc(loci.s3.r.x),...
            loc.s3.y,-paramsCM_opt_nsc(loci.s3.r.z)]';
        locSphere_s3_r_opt = [paramsCM_opt_nsc(loci.s3.r.x),...
            loc.s3.y,paramsCM_opt_nsc(loci.s3.r.z)]';
        locSphere_s4_l_opt = [paramsCM_opt_nsc(loci.s4.r.x),...
            loc.s4.y,-paramsCM_opt_nsc(loci.s4.r.z)]';
        locSphere_s4_r_opt = [paramsCM_opt_nsc(loci.s4.r.x),...
            loc.s4.y,paramsCM_opt_nsc(loci.s4.r.z)]';
        locSphere_s5_l_opt = [paramsCM_opt_nsc(loci.s5.r.x),...
            loc.s5.y,-paramsCM_opt_nsc(loci.s5.r.z)]';
        locSphere_s5_r_opt = [paramsCM_opt_nsc(loci.s5.r.x),...
            loc.s5.y,paramsCM_opt_nsc(loci.s5.r.z)]';
        locSphere_s6_l_opt = [paramsCM_opt_nsc(loci.s6.r.x),...
            loc.s6.y,-paramsCM_opt_nsc(loci.s6.r.z)]';
        locSphere_s6_r_opt = [paramsCM_opt_nsc(loci.s6.r.x),...
            loc.s6.y,paramsCM_opt_nsc(loci.s6.r.z)]';        
        % Compute contact forces
        force_s1_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s1),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s1_l_opt,...
            out1_res_opt(i,calcn.l.pos),out1_res_opt(i,calcn.l.v_lin),...
            out1_res_opt(i,calcn.l.omega),out1_res_opt(i,calcn.TR.R.l),...
            out1_res_opt(i,calcn.TR.T.l));
        force_s2_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s2),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s2_l_opt,...
            out1_res_opt(i,calcn.l.pos),out1_res_opt(i,calcn.l.v_lin),...
            out1_res_opt(i,calcn.l.omega),out1_res_opt(i,calcn.TR.R.l),...
            out1_res_opt(i,calcn.TR.T.l));
        force_s3_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s3),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s3_l_opt,...
            out1_res_opt(i,calcn.l.pos),out1_res_opt(i,calcn.l.v_lin),...
            out1_res_opt(i,calcn.l.omega),out1_res_opt(i,calcn.TR.R.l),...
            out1_res_opt(i,calcn.TR.T.l));
        force_s4_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s4),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s4_l_opt,...
            out1_res_opt(i,toes.l.pos),out1_res_opt(i,toes.l.v_lin),...
            out1_res_opt(i,toes.l.omega),out1_res_opt(i,toes.TR.R.l),...
            out1_res_opt(i,toes.TR.T.l));     
        force_s5_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s5),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s5_l_opt,...
            out1_res_opt(i,calcn.l.pos),out1_res_opt(i,calcn.l.v_lin),...
            out1_res_opt(i,calcn.l.omega),out1_res_opt(i,calcn.TR.R.l),...
            out1_res_opt(i,calcn.TR.T.l));
        force_s6_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s6),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s6_l_opt,...
            out1_res_opt(i,toes.l.pos),out1_res_opt(i,toes.l.v_lin),...
            out1_res_opt(i,toes.l.omega),out1_res_opt(i,toes.TR.R.l),...
            out1_res_opt(i,toes.TR.T.l));        
        force_s1_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s1),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s1_r_opt,...
            out1_res_opt(i,calcn.r.pos),out1_res_opt(i,calcn.r.v_lin),...
            out1_res_opt(i,calcn.r.omega),out1_res_opt(i,calcn.TR.R.r),...
            out1_res_opt(i,calcn.TR.T.r));
        force_s2_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s2),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s2_r_opt,...
            out1_res_opt(i,calcn.r.pos),out1_res_opt(i,calcn.r.v_lin),...
            out1_res_opt(i,calcn.r.omega),out1_res_opt(i,calcn.TR.R.r),...
            out1_res_opt(i,calcn.TR.T.r));
        force_s3_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s3),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s3_r_opt,...
            out1_res_opt(i,calcn.r.pos),out1_res_opt(i,calcn.r.v_lin),...
            out1_res_opt(i,calcn.r.omega),out1_res_opt(i,calcn.TR.R.r),...
            out1_res_opt(i,calcn.TR.T.r));
        force_s4_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s4),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s4_r_opt,...
            out1_res_opt(i,toes.r.pos),out1_res_opt(i,toes.r.v_lin),...
            out1_res_opt(i,toes.r.omega),out1_res_opt(i,toes.TR.R.r),...
            out1_res_opt(i,toes.TR.T.r));   
        force_s5_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s5),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s5_r_opt,...
            out1_res_opt(i,calcn.r.pos),out1_res_opt(i,calcn.r.v_lin),...
            out1_res_opt(i,calcn.r.omega),out1_res_opt(i,calcn.TR.R.r),...
            out1_res_opt(i,calcn.TR.T.r));
        force_s6_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s6),...
            dissipation,normal,transitionVelo,staticFriction,...
            dynamicFriction,viscousFriction,locSphere_s6_r_opt,...
            out1_res_opt(i,toes.r.pos),out1_res_opt(i,toes.r.v_lin),...
            out1_res_opt(i,toes.r.omega),out1_res_opt(i,toes.TR.R.r),...
            out1_res_opt(i,toes.TR.T.r));         
        in_F2_opt = [force_s1_l_opt,force_s2_l_opt,force_s3_l_opt,...
            force_s4_l_opt,force_s5_l_opt,force_s6_l_opt,force_s1_r_opt,...
            force_s2_r_opt,force_s3_r_opt,force_s4_r_opt,force_s5_r_opt,...
            force_s6_r_opt,...
            paramsCM_opt_nsc(loci.s1.r.x),...
            paramsCM_opt_nsc(loci.s1.r.z),...
            paramsCM_opt_nsc(loci.s2.r.x),...
            paramsCM_opt_nsc(loci.s2.r.z),...
            paramsCM_opt_nsc(loci.s3.r.x),...
            paramsCM_opt_nsc(loci.s3.r.z),...
            paramsCM_opt_nsc(loci.s4.r.x),...
            paramsCM_opt_nsc(loci.s4.r.z),...
            paramsCM_opt_nsc(loci.s5.r.x),...
            paramsCM_opt_nsc(loci.s5.r.z),...
            paramsCM_opt_nsc(loci.s6.r.x),...
            paramsCM_opt_nsc(loci.s6.r.z),...
            paramsCM_opt_nsc(radi.s1:radi.s6)'];   
        out2_res = ...
            F2([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)';in_F2_opt']);        
        out2_res_opt(i,:) = full(out2_res);  
    end
    % Optimal joint torques, ground reaction forces and moments
    Tauk_out        = out2_res_opt(:,residualsi);
    GRF_opt_unsc    = out2_res_opt(:,GRFi.all);
    GRM_opt_unsc    = out2_res_opt(:,GRMi.all);
    % assertArmTmax should be 0
    assertArmTmax = max(max(abs(out2_res_opt(:,armsi)-(a_a_opt_unsc)*...
        scaling.ArmTau)));   
   
    %% Passive joint torques at optimal solution
    Tau_pass_opt_all = zeros(N,15);
    for i = 1:N    
        Tau_pass_opt.hip.flex.l    = f_PassiveMoments(k_pass.hip.flex,...
            theta.pass.hip.flex,Xk_Qs_Qdots_opt(i,jointi.hip_flex.l*2-1),...
            Xk_Qs_Qdots_opt(i,jointi.hip_flex.l*2));
        Tau_pass_opt.hip.flex.r    = f_PassiveMoments(k_pass.hip.flex,...
            theta.pass.hip.flex,Xk_Qs_Qdots_opt(i,jointi.hip_flex.r*2-1),...
            Xk_Qs_Qdots_opt(i,jointi.hip_flex.r*2));
        Tau_pass_opt.hip.add.l     = f_PassiveMoments(k_pass.hip.add,...
            theta.pass.hip.add,Xk_Qs_Qdots_opt(i,jointi.hip_add.l*2-1),...
            Xk_Qs_Qdots_opt(i,jointi.hip_add.l*2));
        Tau_pass_opt.hip.add.r     = f_PassiveMoments(k_pass.hip.add,...
            theta.pass.hip.add,Xk_Qs_Qdots_opt(i,jointi.hip_add.r*2-1),...
            Xk_Qs_Qdots_opt(i,jointi.hip_add.r*2));
        Tau_pass_opt.hip.rot.l     = f_PassiveMoments(k_pass.hip.rot,...
            theta.pass.hip.rot,Xk_Qs_Qdots_opt(i,jointi.hip_rot.l*2-1),...
            Xk_Qs_Qdots_opt(i,jointi.hip_rot.l*2));
        Tau_pass_opt.hip.rot.r     = f_PassiveMoments(k_pass.hip.rot,...
            theta.pass.hip.rot,Xk_Qs_Qdots_opt(i,jointi.hip_rot.r*2-1),...
            Xk_Qs_Qdots_opt(i,jointi.hip_rot.r*2));
        Tau_pass_opt.knee.l        = f_PassiveMoments(k_pass.knee,...
            theta.pass.knee,Xk_Qs_Qdots_opt(i,jointi.knee.l*2-1),...
            Xk_Qs_Qdots_opt(i,jointi.knee.l*2));
        Tau_pass_opt.knee.r        = f_PassiveMoments(k_pass.knee,...
            theta.pass.knee,Xk_Qs_Qdots_opt(i,jointi.knee.r*2-1),...
            Xk_Qs_Qdots_opt(i,jointi.knee.r*2));
        Tau_pass_opt.ankle.l       = f_PassiveMoments(k_pass.ankle,...
            theta.pass.ankle,Xk_Qs_Qdots_opt(i,jointi.ankle.l*2-1),...
            Xk_Qs_Qdots_opt(i,jointi.ankle.l*2));
        Tau_pass_opt.ankle.r       = f_PassiveMoments(k_pass.ankle,...
            theta.pass.ankle,Xk_Qs_Qdots_opt(i,jointi.ankle.r*2-1),...
            Xk_Qs_Qdots_opt(i,jointi.ankle.r*2));
        Tau_pass_opt.subt.l       = f_PassiveMoments(k_pass.subt,...
            theta.pass.subt,Xk_Qs_Qdots_opt(i,jointi.subt.l*2-1),...
            Xk_Qs_Qdots_opt(i,jointi.subt.l*2));
        Tau_pass_opt.subt.r       = f_PassiveMoments(k_pass.subt,...
            theta.pass.subt,Xk_Qs_Qdots_opt(i,jointi.subt.r*2-1),...
            Xk_Qs_Qdots_opt(i,jointi.subt.r*2));
        Tau_pass_opt.trunk.ext     = f_PassiveMoments(k_pass.trunk.ext,...
            theta.pass.trunk.ext,Xk_Qs_Qdots_opt(i,jointi.trunk.ext*2-1),...
            Xk_Qs_Qdots_opt(i,jointi.trunk.ext*2));
        Tau_pass_opt.trunk.ben     = f_PassiveMoments(k_pass.trunk.ben,...
            theta.pass.trunk.ben,Xk_Qs_Qdots_opt(i,jointi.trunk.ben*2-1),...
            Xk_Qs_Qdots_opt(i,jointi.trunk.ben*2));
        Tau_pass_opt.trunk.rot     = f_PassiveMoments(k_pass.trunk.rot,...
            theta.pass.trunk.rot,Xk_Qs_Qdots_opt(i,jointi.trunk.rot*2-1),...
            Xk_Qs_Qdots_opt(i,jointi.trunk.rot*2));        
        Tau_pass_opt_all(i,:) = full([Tau_pass_opt.hip.flex.l,...
            Tau_pass_opt.hip.add.l,Tau_pass_opt.hip.rot.l,...             
            Tau_pass_opt.hip.flex.r,Tau_pass_opt.hip.add.r,...
            Tau_pass_opt.hip.rot.r,Tau_pass_opt.knee.l,...
            Tau_pass_opt.knee.r,Tau_pass_opt.ankle.l,...
            Tau_pass_opt.ankle.r,Tau_pass_opt.subt.l,...
            Tau_pass_opt.subt.r,Tau_pass_opt.trunk.ext,...
            Tau_pass_opt.trunk.ben,Tau_pass_opt.trunk.rot]);
    end        
    
    %% Create .mot file for OpenSim GUI
    if writeIKmotion  
        pathOpenSim = [pathRepo,'/OpenSim'];
        addpath(genpath(pathOpenSim));
        q_opt_GUI = zeros(N+1,1+nq.all+2);
        q_opt_GUI(:,1) = tgrid';
        q_opt_GUI(:,2:nq.all+1)  = q_opt_unsc_all.deg;
        q_opt_GUI(:,end-1:end) = 1.51*180/pi*ones(N+1,2);% pro_sup (locked)
        JointAngle.labels = {'time','pelvis_tilt','pelvis_list',...
            'pelvis_rotation','pelvis_tx','pelvis_ty','pelvis_tz',...
            'hip_flexion_l','hip_adduction_l','hip_rotation_l',...
            'hip_flexion_r','hip_adduction_r','hip_rotation_r',...
            'knee_angle_l','knee_angle_r','ankle_angle_l',...
            'ankle_angle_r','subtalar_angle_l','subtalar_angle_r',...
            'lumbar_extension','lumbar_bending','lumbar_rotation',...
            'arm_flex_l','arm_add_l','arm_rot_l',...
            'arm_flex_r','arm_add_r','arm_rot_r',...
            'elbow_flex_l','elbow_flex_r',...
            'pro_sup_l','pro_sup_r'};
        % Combine data joint angles and muscle activations
        JointAngleMuscleAct.data = [q_opt_GUI,a_opt];
        % Get muscle labels
        muscleNamesAll = cell(1,NMuscle);
        for i = 1:NMuscle/2
            muscleNamesAll{i} = [muscleNames{i}(1:end-2),'_l'];
            muscleNamesAll{i+NMuscle/2} = [muscleNames{i}(1:end-2),'_r'];
        end  
        % Combine labels joint angles and muscle activations
        JointAngleMuscleAct.labels = JointAngle.labels;
        for i = 1:NMuscle
            JointAngleMuscleAct.labels{i+size(q_opt_GUI,2)} = ...
                [muscleNamesAll{i},'/activation'];
        end        
        filenameJointAngles = [pathRepo,'/Results/',namescript,...
                '/IK',savename,'.mot'];
        write_motionFile(JointAngleMuscleAct, filenameJointAngles)
    end   
    
    %% Optimal cost and CPU time
    pathDiary = [pathresults,'/',namescript,'/D',savename];
    [CPU_IPOPT,CPU_NLP,~,~,~,~,~,~,OptSol] = readDiary(pathDiary);
        
    %% Save results       
    if saveResults
        if (exist([pathresults,'/',namescript,'/Results_tracking.mat'],...
            'file')==2) 
        load([pathresults,'/',namescript,'/Results_tracking.mat']);
    else
        Results_tracking = struct('Qs_opt',[]);
        end
    % Structure results
    Results_tracking.Qs_opt = q_opt_unsc.deg;
    Results_tracking.Acts_opt = a_opt_unsc;
    Results_tracking.Ts_opt = Tauk_out;
    Results_tracking.GRFs_opt = GRF_opt_unsc;
    Results_tracking.GRMs_opt = GRM_opt_unsc; 
    Results_tracking.ParamsCM_opt = full(paramsCM_opt_nsc);     
    Results_tracking.Qs_toTrack = Qs.allinterpfilt;
    Results_tracking.Ts_toTrack = ID.allinterp;
    Results_tracking.GRFs_toTrack = GRF.val.allinterp;
    Results_tracking.GRMs_toTrack = GRF.MorGF.allinterp; 
    Results_tracking.ParamsCM_gen = ...
        full(f_nsc18(guess.params,scaling.params.v,scaling.params.r));
    Results_tracking.CPU_IPOPT = CPU_IPOPT;
    Results_tracking.CPU_NLP = CPU_NLP;
    Results_tracking.OptSol = OptSol;
    Results_tracking.colheaders.joints = joints;
    Results_tracking.colheaders.GRF = {'fore_aft_r','vertical_r',...
        'lateral_r','fore_aft_l','vertical_l','lateral_l'};
    for i = 1:NMuscle/2
            Results_tracking.colheaders.muscles{i} = ...
                [muscleNames{i}(1:end-2),'_l'];
            Results_tracking.colheaders.muscles{i+NMuscle/2} = ...
                [muscleNames{i}(1:end-2),'_r'];
    end
    Results_tracking.colheaders.paramsCM = {'loc_s1_r_x','loc_s1_r_z',...
        'loc_s2_r_x','loc_s2_r_z','loc_s3_r_x','loc_s3_r_z',...
        'loc_s4_r_x','loc_s4_r_z','loc_s5_r_x','loc_s5_r_z',...
        'loc_s6_r_x','loc_s6_r_z','radius_s1','radius_s2','radius_s3',...
        'radius_s4','radius_s5','radius_s6'};
    % Save data
    save([pathresults,'/',namescript,'/Results_tracking.mat'],...
        'Results_tracking');
    end
end
end
