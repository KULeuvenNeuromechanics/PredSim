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
% num_set = [0,1,1,1,0,0]; % This configuration analyzes the results

% The variable settings in the following section will set some parameters 
% of the optimization problem. Through the variable idx_ww, the user can 
% select which row of parameters will be used.
idx_ww = 25; % Index row in matrix settings

%% Settings
import casadi.*
subjectData = 'subject1';

solveProblem    = num_set(1); % set to 1 to solve problem
analyseResults  = num_set(2); % set to 1 to analyze results
loadResults     = num_set(3); % set to 1 to load results
saveResults     = num_set(4); % set to 1 to save sens. results
checkBoundsIG   = num_set(5); % set to 1 to visualize guess-bounds 
writeIKmotion   = num_set(6); % set to 1 to write .mot file

% settings describes the parameters used in the optimal control problem.
trackSim_mtp_settings_all

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
cs          = settings(ww,10); % number of contact spheres
csc         = settings(ww,11); % contact sphere configuration
sd          = settings(ww,12);  % stiffness and damping
% Trials
trials = settings_trials(ww).ww;
% Fixed parameter
W.u = 0.001;
% Identifiers for experimental data

if (cs == 5 && csc == 3) || (cs == 4 && csc == 3) || (cs == 6 && csc == 4) || (cs == 6 && csc == 5) || (cs == 6 && csc == 6)
    mtp_jointType = 'mtpPin';
else
    mtp_jointType = 'mtp';
end
subject = ['subject1_',mtp_jointType];

for p = 1:length(trials)
    nametrial(p).id    = ['gait_',trials{p},'_',mtp_jointType]; % Experimental walking trial to track
    nametrial(p).ID    = ['ID_',nametrial(p).id];
    nametrial(p).GRF   = ['GRF_',nametrial(p).id];
    nametrial(p).IK    = ['IK_',nametrial(p).id];
    switch nametrial(p).id
        case ['gait_14_',mtp_jointType]
            time_opt(p).p = [3.73,4.25];
        case ['gait_15_',mtp_jointType]
            time_opt(p).p = [2.66,3.2];
    end  
end
id_case = strjoin(trials,'-');
% The filename used to save the results depends on the settings 
savename = ['_c',num2str(ww)];

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
            if (cs == 5 && csc == 3) || (cs == 4 && csc == 3) || (cs == 6 && csc == 4) || (cs == 6 && csc == 5) || (cs == 6 && csc == 6)
                F1 = external('F','TrackSim_mtpPin_1.dll');
            else
                F1 = external('F','TrackSim_mtp_1.dll');
            end
            if cs == 6
                if csc == 4 || csc == 5 || csc == 6
                    F2 = external('F','TrackSim_mtpPin_2_6cs.dll'); 
                else
                    F2 = external('F','TrackSim_mtp_2.dll'); 
                end
            elseif cs == 5
                if csc == 1
                    F2 = external('F','TrackSim_mtp_2_5cs.dll');  
                elseif csc == 2
                    F2 = external('F','TrackSim_mtp_2_5csb.dll'); 
                elseif csc == 3
                    F2 = external('F','TrackSim_mtpPin_2_5csb.dll');
                end
            elseif cs == 4
                if csc == 1
                    F2 = external('F','TrackSim_mtp_2_4cs.dll'); 
                elseif csc == 2
                    F2 = external('F','TrackSim_mtp_2_4csb.dll'); 
                elseif csc == 3
                    F2 = external('F','TrackSim_mtpPin_2_4cs.dll');
                end
            end
    end
elseif ismac
    switch setup.derivatives
        case 'AD'     
            cd(pathExternalFunctions);
            F1 = external('F','TrackSim_1.dylib');
            if cs == 6                 
                F2 = external('F','TrackSim_2.dylib'); 
            else
                disp('Platform not supported')
                continue
            end
    end
else
    disp('Platform not supported')
    continue
end   
cd(pathmain);
% This is an example of how to call an external function with some
% numerical values.
% vec1 = -ones(58,1);
% res1 = full(F3(vec1));
% vec2 = -ones(62,1);
% vec2([19*2-1:20*2]) = 0;
% res2 = full(F1(vec2));
% max(abs(res1)-abs(res2))
% vec1 = -ones(141,1);
% res1 = full(F4(vec1));
% vec2 = -ones(147,1);
% vec2([19*2-1:20*2, 31*2+19:31*2+20]) = 0;
% res2 = full(F2(vec2));
% max(abs(res1)-abs(res2([1:18,21:end])))

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
mtpi                = jointi.mtp.l:jointi.mtp.r; % mtp
residuals_acti      = jointi.hip_flex.l:jointi.elb.r; % all but gr-pelvis
residuals_act_bmtpi = [jointi.hip_flex.l:jointi.subt.r,...
    jointi.trunk.ext:jointi.elb.r];
residual_bptyi      = [jointi.pelvis.tilt:jointi.pelvis.tx,...
    jointi.pelvis.tz:jointi.elb.r]; % all but pelvis_ty
% Number of degrees of freedom for later use
nq.all              = length(residualsi); % all 
nq.abs              = length(ground_pelvisi); % ground-pelvis
nq.act              = nq.all-nq.abs;% all but ground-pelvis
nq.trunk            = length(trunki); % trunk
nq.arms             = length(armsi); % arms
nq.mtp              = length(mtpi); % mtp
nq.leg              = 10; % #joints needed for polynomials
Qsi                 = 1:2:2*nq.all; % indices Qs only
% Second, GRFs
GRFi.r              = 32:34;
GRFi.l              = 35:37;
GRFi.all            = [GRFi.r,GRFi.l];
nGRF                = length(GRFi.all);
% Third, GRMs
GRMi.r              = 38:40;
GRMi.l              = 41:43;
GRMi.all            = [GRMi.r,GRMi.l];
nGRM                = length(GRMi.all);
% Number contact model parameters
np = cs*3; % transversal location (2) + radius (1) per contact sphere  

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
% Fixed parameters
dissipation = 2;
normal = [0,1,0];
transitionVelo = 0.2;
staticFriction = 0.8;
dynamicFriction = 0.8;
viscousFriction = 0.5; 
stif = 1000000;  
if (cs == 6 && csc == 4) || (cs == 6 && csc == 5) || (cs == 6 && csc == 6)
    loc.s1.y = -0.01;
    loc.s2.y = -0.01;
    loc.s3.y = -0.01;
    loc.s4.y = -0.01;
    loc.s5.y = -0.01;
    loc.s6.y = -0.01;
else
    loc.s1.y = -0.021859;
    loc.s2.y = -0.021859;
    loc.s3.y = -0.021859;
    loc.s4.y = -0.0214476;
    if (cs == 4 && csc == 2) || (cs == 4 && csc == 3) % Sphere 3 => toes
        loc.s3.y = -0.0214476;    
    end
    if cs == 5 
        if csc == 1 % Sphere 5 => calcn
            loc.s5.y = -0.021859;    
        elseif csc == 2 || csc == 3 % Sphere 5 => toes
            loc.s5.y = -0.0214476;   
        end
    end
    if cs == 6 
        loc.s6.y = -0.0214476;
    end
end
% Indices variable parameters
loci.s1.r.x = 1;
loci.s1.r.z = 2;
loci.s2.r.x = 3;
loci.s2.r.z = 4;
loci.s3.r.x = 5;
loci.s3.r.z = 6;
loci.s4.r.x = 7;
loci.s4.r.z = 8;
radi.s1 = 2*cs+1;
radi.s2 = 2*cs+2;
radi.s3 = 2*cs+3;
radi.s4 = 2*cs+4;
if cs == 5    
    loci.s5.r.x = 2*cs-1;
    loci.s5.r.z = 2*cs;
    radi.s5 = 2*cs+5;
elseif cs == 6        
    loci.s5.r.x = 2*cs-3;
    loci.s5.r.z = 2*cs-2;
    loci.s6.r.x = 2*cs-1;
    loci.s6.r.z = 2*cs;    
    radi.s5 = 2*cs+5;
    radi.s6 = 2*cs+6;
end

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
CasADiFunctions_tracking_mtp

%% Passive joint torques
% We extract the parameters for the passive torques of the lower limbs and
% the trunk
pathPassiveMoments = [pathRepo,'/PassiveMoments'];
addpath(genpath(pathPassiveMoments));
PassiveMomentsData
stiffnessArm = 0;
dampingArm = 0.1;
stiffnessMtp = 1.5/(pi/180)/5;
dampingMtp = 0.5/5;
if sd == 2
    stiffnessMtp = 1.5/(pi/180);
    dampingMtp = 0.5;
elseif sd == 3
    stiffnessMtp = 1.5/(pi/180)/2;
    dampingMtp = 0.5/2;
elseif sd == 4
    stiffnessMtp = 1.5/(pi/180)/5;
    dampingMtp = 0.5;
elseif sd == 5
    stiffnessMtp = 1.5/(pi/180)/3;
    dampingMtp = 0.5;
elseif sd == 6
    stiffnessMtp = 1.5/(pi/180)/5;
    dampingMtp = 0.5;
end    

%% Experimental data
pathData = [pathRepo,'/OpenSimModel/',subjectData];
joints = {'pelvis_tilt','pelvis_list','pelvis_rotation','pelvis_tx',...
    'pelvis_ty','pelvis_tz','hip_flexion_l','hip_adduction_l',...
    'hip_rotation_l','hip_flexion_r','hip_adduction_r','hip_rotation_r',...
    'knee_angle_l','knee_angle_r','ankle_angle_l','ankle_angle_r',...
    'subtalar_angle_l','subtalar_angle_r',...
    'mtp_angle_l','mtp_angle_r', ...
    'lumbar_extension','lumbar_bending','lumbar_rotation',...
    'arm_flex_l','arm_add_l','arm_rot_l',...
    'arm_flex_r','arm_add_r','arm_rot_r',...
    'elbow_flex_l','elbow_flex_r'};
pathVariousFunctions = [pathRepo,'/VariousFunctions'];
addpath(genpath(pathVariousFunctions));

for p = 1:length(trials)
    % Extract joint kinematics
    pathIK = [pathData,'/IK/',nametrial(p).IK,'.mat'];
    Qs(p).p = getIK(pathIK,joints);
    % Extract ground reaction forces and moments
    pathGRF = [pathData,'/GRF/',nametrial(p).GRF,'.mat'];
    GRF(p).p = getGRF(pathGRF);
    % Extract joint kinetics
    pathID = [pathData,'/ID/',nametrial(p).ID,'.mat'];
    ID(p).p = getID(pathID,joints);
    % Interpolation experimental data
    time_expi.ID(1) = find(round(ID(p).p.time,4) == time_opt(p).p (1));
    time_expi.ID(2) = find(round(ID(p).p.time,4) == time_opt(p).p (2));
    time_expi.GRF(1) = find(round(GRF(p).p.time,4) == time_opt(p).p (1));
    time_expi.GRF(2) = find(round(GRF(p).p.time,4) == time_opt(p).p (2));
    step = (ID(p).p.time(time_expi.ID(2))-ID(p).p.time(time_expi.ID(1)))/(N-1);
    interval = time_opt(p).p (1):step:time_opt(p).p (2);
    ID(p).p.allinterp = interp1(ID(p).p.all(:,1),ID(p).p.all,interval);
    Qs(p).p.allinterpfilt = interp1(Qs(p).p.allfilt(:,1),Qs(p).p.allfilt,interval);
    GRF(p).p.val.allinterp = interp1(round(GRF(p).p.val.all(:,1),4),...
        GRF(p).p.val.all,round(interval,4));
    GRF(p).p.MorGF.allinterp = interp1(round(GRF(p).p.MorGF.all(:,1),4),...
        GRF(p).p.MorGF.all,round(interval,4));
    GRF(p).p.pos.allinterp = interp1(round(GRF(p).p.pos.all(:,1),4),...
        GRF(p).p.pos.all,round(interval,4));
    GRF(p).p.Mcop.allinterp = interp1(round(GRF(p).p.Mcop.all(:,1),4),...
        GRF(p).p.Mcop.all,round(interval,4));
end

%% Bounds
pathBounds = [pathRepo,'/Bounds'];
addpath(genpath(pathBounds));
for p = 1:length(trials)
    [bounds(p).p,scaling(p).p] = ...
        getBounds_tracking_mtp(Qs(p).p,NMuscle,nq,jointi,dev_cm,GRF(p).p,cs,csc);
end

%% Initial guess
pathIG = [pathRepo,'/IG'];
addpath(genpath(pathIG));
% Data-informed initial guess
for p = 1:length(trials)
    guess(p).p = ...
        getGuess_DI_tracking_mtp(Qs(p).p,nq,N,NMuscle,jointi,scaling(p).p,cs,csc);
end
% This allows visualizing the initial guess and the boundsbound
if checkBoundsIG
    pathPlots = [pathRepo,'/Plots'];
    addpath(genpath(pathPlots));
    for p = 1:length(trials)
        plot_BoundsVSInitialGuess_tracking_mtp_mp
    end
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
    paramsCM     = MX.sym('paramsCM',np);
    w            = [w {paramsCM}];
    lbw          = [lbw; bounds(1).p.params.lower'];
    ubw          = [ubw; bounds(1).p.params.upper'];
    w0           = [w0;  guess(1).p.params'];      
    paramsCM_nsc = f_nscnp(paramsCM,scaling(1).p.params.v,scaling(1).p.params.r);        
    for p = 1:length(trials)    
        % Define states at first mesh point
        % Muscle activations
        a0              = MX.sym('a0',NMuscle);
        w               = [w {a0}];
        lbw             = [lbw; bounds(p).p.a.lower'];
        ubw             = [ubw; bounds(p).p.a.upper'];
        w0              = [w0;  guess(p).p.a(1,:)'];
        % Muscle-tendon forces
        FTtilde0        = MX.sym('FTtilde0',NMuscle);
        w               = [w {FTtilde0}];
        lbw             = [lbw; bounds(p).p.FTtilde.lower'];
        ubw             = [ubw; bounds(p).p.FTtilde.upper'];
        w0              = [w0;  guess(p).p.FTtilde(1,:)'];
        % Qs and Qdots
        X0              = MX.sym('X0',2*nq.all);
        w               = [w {X0}];    
        lbw             = [lbw; bounds(p).p.QsQdots.lower'];
        ubw             = [ubw; bounds(p).p.QsQdots.upper'];    
        w0              = [w0;  guess(p).p.QsQdots(1,:)'];
        % Arm activations
        a_a0            = MX.sym('a_a0',nq.arms);
        w               = [w {a_a0}];
        lbw             = [lbw; bounds(p).p.a_a.lower'];
        ubw             = [ubw; bounds(p).p.a_a.upper'];
        w0              = [w0;  guess(p).p.a_a(1,:)'];  
        % Mtp activations
        a_mtp0          = MX.sym('a_mtp0',nq.mtp);
        w               = [w {a_mtp0}];
        lbw             = [lbw; bounds(p).p.a_mtp.lower'];
        ubw             = [ubw; bounds(p).p.a_mtp.upper'];
        w0              = [w0;  guess(p).p.a_mtp(1,:)'];  
        % "Lift" initial conditions
        ak          = a0;
        FTtildek    = FTtilde0;
        Xk          = X0;
        a_ak        = a_a0;
        a_mtpk      = a_mtp0;
        % Time step
        h = (time_opt(p).p (2)-time_opt(p).p (1))/N;
        % Loop over mesh points
        for k=0:N-1
            % Define controls at mesh point (piecewise-constant in interval) 
            % Time derivative of muscle activations (states)
            vAk                 = MX.sym(['vA_' num2str(k)], NMuscle);
            w                   = [w {vAk}];
            lbw                 = [lbw; bounds(p).p.vA.lower'];
            ubw                 = [ubw; bounds(p).p.vA.upper'];
            w0                  = [w0; guess(p).p.vA(k+1,:)'];
            % Time derivative of muscle-tendon forces (states)
            dFTtildek           = MX.sym(['dFTtilde_' num2str(k)], NMuscle);
            w                   = [w {dFTtildek}];
            lbw                 = [lbw; bounds(p).p.dFTtilde.lower'];
            ubw                 = [ubw; bounds(p).p.dFTtilde.upper'];
            w0                  = [w0; guess(p).p.dFTtilde(k+1,:)'];  
            % Time derivative of Qdots (states) 
            Ak                  = MX.sym(['A_' num2str(k)], nq.all);
            w                   = [w {Ak}];
            lbw                 = [lbw; bounds(p).p.Qdotdots.lower'];
            ubw                 = [ubw; bounds(p).p.Qdotdots.upper'];
            w0                  = [w0; guess(p).p.Qdotdots(k+1,:)'];    
            % Arm excitations
            e_ak                = MX.sym(['e_a_' num2str(k)], nq.arms);
            w                   = [w {e_ak}];
            lbw                 = [lbw; bounds(p).p.e_a.lower'];
            ubw                 = [ubw; bounds(p).p.e_a.upper'];
            w0                  = [w0; guess(p).p.e_a(k+1,:)'];
            % Mtp excitations
            e_mtpk              = MX.sym(['e_mtp_' num2str(k)], nq.mtp);
            w                   = [w {e_mtpk}];
            lbw                 = [lbw; bounds(p).p.e_mtp.lower'];
            ubw                 = [ubw; bounds(p).p.e_mtp.upper'];
            w0                  = [w0; guess(p).p.e_mtp(k+1,:)'];
            % Define states at collocation points    
            % Muscle activations
            akj = {};
            for j=1:d
                akj{j}  = MX.sym(['	a_' num2str(k) '_' num2str(j)], NMuscle);
                w       = {w{:}, akj{j}};
                lbw     = [lbw; bounds(p).p.a.lower'];
                ubw     = [ubw; bounds(p).p.a.upper'];
                w0      = [w0;  guess(p).p.a(k+1,:)'];
            end   
            % Muscle-tendon forces
            FTtildekj = {};
            for j=1:d
                FTtildekj{j} = ...
                    MX.sym(['FTtilde_' num2str(k) '_' num2str(j)], NMuscle);
                w            = {w{:}, FTtildekj{j}};
                lbw          = [lbw; bounds(p).p.FTtilde.lower'];
                ubw          = [ubw; bounds(p).p.FTtilde.upper'];
                w0           = [w0;  guess(p).p.FTtilde(k+1,:)'];
            end
            % Qs and Qdots        
            Xkj = {};
            for j=1:d
                Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], 2*nq.all);
                w      = {w{:}, Xkj{j}};
                lbw    = [lbw; bounds(p).p.QsQdots.lower'];
                ubw    = [ubw; bounds(p).p.QsQdots.upper'];
                w0     = [w0;  guess(p).p.QsQdots(k+1,:)'];
            end   
            % Arm activations
            a_akj = {};
            for j=1:d
                a_akj{j}= MX.sym(['a_a_' num2str(k) '_' num2str(j)], nq.arms);
                w       = {w{:}, a_akj{j}};
                lbw     = [lbw; bounds(p).p.a_a.lower'];
                ubw     = [ubw; bounds(p).p.a_a.upper'];
                w0      = [w0;  guess(p).p.a_a(k+1,:)'];
            end   
            % Mtp activations
            a_mtpkj = {};
            for j=1:d
                a_mtpkj{j}= MX.sym(['a_mtp_' num2str(k) '_' num2str(j)], nq.mtp);
                w       = {w{:}, a_mtpkj{j}};
                lbw     = [lbw; bounds(p).p.a_mtp.lower'];
                ubw     = [ubw; bounds(p).p.a_mtp.upper'];
                w0      = [w0;  guess(p).p.a_mtp(k+1,:)'];
            end   
            % Unscale variables for later use
            Xk_nsc          = Xk.*scaling(p).p.QsQdots';
            FTtildek_nsc    = FTtildek.*(scaling(p).p.FTtilde');
            Ak_nsc          = Ak.*scaling(p).p.Qdotdots';
            for j=1:d
                Xkj_nsc{j} = Xkj{j}.*scaling(p).p.QsQdots';
                FTtildekj_nsc{j} = FTtildekj{j}.*scaling(p).p.FTtilde';
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
            if cs == 5 || cs == 6
                locSphere_s5_l = [paramsCM_nsc(loci.s5.r.x),...
                    loc.s5.y,-paramsCM_nsc(loci.s5.r.z)]';
                locSphere_s5_r = [paramsCM_nsc(loci.s5.r.x),...
                    loc.s5.y,paramsCM_nsc(loci.s5.r.z)]';     
                if cs == 6
                    locSphere_s6_l = [paramsCM_nsc(loci.s6.r.x),...
                        loc.s6.y,-paramsCM_nsc(loci.s6.r.z)]';
                    locSphere_s6_r = [paramsCM_nsc(loci.s6.r.x),...
                        loc.s6.y,paramsCM_nsc(loci.s6.r.z)]'; 
                end
            end
            % Compute contact forces            
            force_s1_l = f_contactForce(stif,paramsCM_nsc(radi.s1),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s1_l,...
                outk1(calcn.l.pos),outk1(calcn.l.v_lin),...
                outk1(calcn.l.omega),outk1(calcn.TR.R.l),outk1(calcn.TR.T.l));
            force_s1_r = f_contactForce(stif,paramsCM_nsc(radi.s1),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s1_r,...
                outk1(calcn.r.pos),outk1(calcn.r.v_lin),...
                outk1(calcn.r.omega),outk1(calcn.TR.R.r),outk1(calcn.TR.T.r));
            force_s2_l = f_contactForce(stif,paramsCM_nsc(radi.s2),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s2_l,...
                outk1(calcn.l.pos),outk1(calcn.l.v_lin),...
                outk1(calcn.l.omega),outk1(calcn.TR.R.l),outk1(calcn.TR.T.l));                     
            force_s2_r = f_contactForce(stif,paramsCM_nsc(radi.s2),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s2_r,...
                outk1(calcn.r.pos),outk1(calcn.r.v_lin),...
                outk1(calcn.r.omega),outk1(calcn.TR.R.r),outk1(calcn.TR.T.r)); 
            if (cs == 4 && csc == 2) || (cs == 4 && csc == 3) % Sphere 3 => toes            
                force_s3_l = f_contactForce(stif,paramsCM_nsc(radi.s3),...
                    dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,locSphere_s3_l,...
                    outk1(toes.l.pos),outk1(toes.l.v_lin),...
                    outk1(toes.l.omega),outk1(toes.TR.R.l),outk1(toes.TR.T.l));         
                force_s3_r = f_contactForce(stif,paramsCM_nsc(radi.s3),...
                    dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,locSphere_s3_r,...
                    outk1(toes.r.pos),outk1(toes.r.v_lin),...
                    outk1(toes.r.omega),outk1(toes.TR.R.r),outk1(toes.TR.T.r));     
            else % Sphere 3 => calc 
                force_s3_l = f_contactForce(stif,paramsCM_nsc(radi.s3),...
                    dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,locSphere_s3_l,...
                    outk1(calcn.l.pos),outk1(calcn.l.v_lin),...
                    outk1(calcn.l.omega),outk1(calcn.TR.R.l),outk1(calcn.TR.T.l));
                force_s3_r = f_contactForce(stif,paramsCM_nsc(radi.s3),...
                    dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,locSphere_s3_r,...
                    outk1(calcn.r.pos),outk1(calcn.r.v_lin),...
                    outk1(calcn.r.omega),outk1(calcn.TR.R.r),outk1(calcn.TR.T.r));
            end
            if (cs == 6 && csc == 4) || (cs == 6 && csc == 5) || (cs == 6 && csc == 6) % Sphere 4 => calcn 
                force_s4_l = f_contactForce(stif,paramsCM_nsc(radi.s4),...
                    dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,locSphere_s4_l,...
                    outk1(calcn.l.pos),outk1(calcn.l.v_lin),...
                    outk1(calcn.l.omega),outk1(calcn.TR.R.l),outk1(calcn.TR.T.l)); 
                force_s4_r = f_contactForce(stif,paramsCM_nsc(radi.s4),...
                    dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,locSphere_s4_r,...
                    outk1(calcn.r.pos),outk1(calcn.r.v_lin),...
                    outk1(calcn.r.omega),outk1(calcn.TR.R.r),outk1(calcn.TR.T.r)); 
            else  % Sphere 4 => toes 
                force_s4_l = f_contactForce(stif,paramsCM_nsc(radi.s4),...
                    dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,locSphere_s4_l,...
                    outk1(toes.l.pos),outk1(toes.l.v_lin),...
                    outk1(toes.l.omega),outk1(toes.TR.R.l),outk1(toes.TR.T.l)); 
                force_s4_r = f_contactForce(stif,paramsCM_nsc(radi.s4),...
                    dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,locSphere_s4_r,...
                    outk1(toes.r.pos),outk1(toes.r.v_lin),...
                    outk1(toes.r.omega),outk1(toes.TR.R.r),outk1(toes.TR.T.r)); 
            end           
            
            if cs == 5 || cs == 6
                if csc == 1 % Sphere 5 => calcn
                    force_s5_l = f_contactForce(stif,paramsCM_nsc(radi.s5),...
                        dissipation,normal,transitionVelo,staticFriction,...
                        dynamicFriction,viscousFriction,locSphere_s5_l,...
                        outk1(calcn.l.pos),outk1(calcn.l.v_lin),...
                        outk1(calcn.l.omega),outk1(calcn.TR.R.l),...
                        outk1(calcn.TR.T.l)); 
                    force_s5_r = f_contactForce(stif,paramsCM_nsc(radi.s5),...
                        dissipation,normal,transitionVelo,staticFriction,...
                        dynamicFriction,viscousFriction,locSphere_s5_r,...
                        outk1(calcn.r.pos),outk1(calcn.r.v_lin),...
                        outk1(calcn.r.omega),outk1(calcn.TR.R.r),...
                        outk1(calcn.TR.T.r)); 
                elseif csc == 2 || csc == 3 || csc == 4 || csc == 5 || csc == 6 % Sphere 5 => toes
                    force_s5_l = f_contactForce(stif,paramsCM_nsc(radi.s5),...
                        dissipation,normal,transitionVelo,staticFriction,...
                        dynamicFriction,viscousFriction,locSphere_s5_l,...
                        outk1(toes.l.pos),outk1(toes.l.v_lin),...
                        outk1(toes.l.omega),outk1(toes.TR.R.l),...
                        outk1(toes.TR.T.l)); 
                    force_s5_r = f_contactForce(stif,paramsCM_nsc(radi.s5),...
                        dissipation,normal,transitionVelo,staticFriction,...
                        dynamicFriction,viscousFriction,locSphere_s5_r,...
                        outk1(toes.r.pos),outk1(toes.r.v_lin),...
                        outk1(toes.r.omega),outk1(toes.TR.R.r),...
                        outk1(toes.TR.T.r));
                end
                if cs == 6
                    force_s6_l = f_contactForce(stif,paramsCM_nsc(radi.s6),...
                        dissipation,normal,transitionVelo,staticFriction,...
                        dynamicFriction,viscousFriction,locSphere_s6_l,...
                        outk1(toes.l.pos),outk1(toes.l.v_lin),...
                        outk1(toes.l.omega),outk1(toes.TR.R.l),...
                        outk1(toes.TR.T.l));  
                    force_s6_r = f_contactForce(stif,paramsCM_nsc(radi.s6),...
                        dissipation,normal,transitionVelo,staticFriction,...
                        dynamicFriction,viscousFriction,locSphere_s6_r,...
                        outk1(toes.r.pos),outk1(toes.r.v_lin),...
                        outk1(toes.r.omega),outk1(toes.TR.R.r),...
                        outk1(toes.TR.T.r));
                end
            end   
            if cs == 4
                in_F2 = [force_s1_l,force_s2_l,force_s3_l,force_s4_l,...
                    force_s1_r,force_s2_r,force_s3_r,force_s4_r,...       
                    paramsCM_nsc(loci.s1.r.x),paramsCM_nsc(loci.s1.r.z),...
                    paramsCM_nsc(loci.s2.r.x),paramsCM_nsc(loci.s2.r.z),...
                    paramsCM_nsc(loci.s3.r.x),paramsCM_nsc(loci.s3.r.z),...
                    paramsCM_nsc(loci.s4.r.x),paramsCM_nsc(loci.s4.r.z),...
                    paramsCM_nsc(radi.s1:radi.s4)'];
            elseif cs == 5
                in_F2 = [force_s1_l,force_s2_l,force_s3_l,force_s4_l,...
                    force_s5_l,force_s1_r,force_s2_r,force_s3_r,...
                    force_s4_r,force_s5_r,...       
                    paramsCM_nsc(loci.s1.r.x),paramsCM_nsc(loci.s1.r.z),...
                    paramsCM_nsc(loci.s2.r.x),paramsCM_nsc(loci.s2.r.z),...
                    paramsCM_nsc(loci.s3.r.x),paramsCM_nsc(loci.s3.r.z),...
                    paramsCM_nsc(loci.s4.r.x),paramsCM_nsc(loci.s4.r.z),...
                    paramsCM_nsc(loci.s5.r.x),paramsCM_nsc(loci.s5.r.z),...
                    paramsCM_nsc(radi.s1:radi.s5)'];
            elseif cs == 6
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
            end
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
                Xk_nsc(jointi.mtp.l*2-1,1),...
                Xk_nsc(jointi.trunk.ext*2-1,1),...
                Xk_nsc(jointi.trunk.ben*2-1,1),...
                Xk_nsc(jointi.trunk.rot*2-1,1)];  
            qdotin_l = [Xk_nsc(jointi.hip_flex.l*2,1),...
                Xk_nsc(jointi.hip_add.l*2,1),...
                Xk_nsc(jointi.hip_rot.l*2,1),...
                Xk_nsc(jointi.knee.l*2,1),...
                Xk_nsc(jointi.ankle.l*2,1),...
                Xk_nsc(jointi.subt.l*2,1),...
                Xk_nsc(jointi.mtp.l*2,1),...
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
            MA.trunk_ext    =  MA_l([47:49,mai(8).mus.l]',8);
            MA.trunk_ben    =  MA_l([47:49,mai(9).mus.l]',9);
            MA.trunk_rot    =  MA_l([47:49,mai(10).mus.l]',10);
            % Right leg
            qin_r = [Xk_nsc(jointi.hip_flex.r*2-1,1),...
                Xk_nsc(jointi.hip_add.r*2-1,1),...
                Xk_nsc(jointi.hip_rot.r*2-1,1),...
                Xk_nsc(jointi.knee.r*2-1,1),...
                Xk_nsc(jointi.ankle.r*2-1,1),...
                Xk_nsc(jointi.subt.r*2-1,1),...
                Xk_nsc(jointi.mtp.r*2-1,1),...
                Xk_nsc(jointi.trunk.ext*2-1,1),...
                Xk_nsc(jointi.trunk.ben*2-1,1),...
                Xk_nsc(jointi.trunk.rot*2-1,1)];  
            qdotin_r = [Xk_nsc(jointi.hip_flex.r*2,1),...
                Xk_nsc(jointi.hip_add.r*2,1),...
                Xk_nsc(jointi.hip_rot.r*2,1),...
                Xk_nsc(jointi.knee.r*2,1),...
                Xk_nsc(jointi.ankle.r*2,1),...
                Xk_nsc(jointi.subt.r*2,1),...
                Xk_nsc(jointi.mtp.r*2,1),...
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
                    ak,FTtildek.*scaling(p).p.FTtilde',...
                    dFTtildek.*scaling(p).p.dFTtilde,lMTk_lr,vMTk_lr,tensions);
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

            Tau_passk.sh_flex.l = f_passiveTATorques(stiffnessArm, dampingArm, ...
                Xk_nsc(jointi.sh_flex.l*2-1,1), Xk_nsc(jointi.sh_flex.l*2,1));
            Tau_passk.sh_add.l = f_passiveTATorques(stiffnessArm, dampingArm, ...
                Xk_nsc(jointi.sh_add.l*2-1,1), Xk_nsc(jointi.sh_add.l*2,1));
            Tau_passk.sh_rot.l = f_passiveTATorques(stiffnessArm, dampingArm, ...
                Xk_nsc(jointi.sh_rot.l*2-1,1), Xk_nsc(jointi.sh_rot.l*2,1));
            Tau_passk.sh_flex.r = f_passiveTATorques(stiffnessArm, dampingArm, ...
                Xk_nsc(jointi.sh_flex.r*2-1,1), Xk_nsc(jointi.sh_flex.r*2,1));
            Tau_passk.sh_add.r = f_passiveTATorques(stiffnessArm, dampingArm, ...
                Xk_nsc(jointi.sh_add.r*2-1,1), Xk_nsc(jointi.sh_add.r*2,1));
            Tau_passk.sh_rot.r = f_passiveTATorques(stiffnessArm, dampingArm, ...
                Xk_nsc(jointi.sh_rot.r*2-1,1), Xk_nsc(jointi.sh_rot.r*2,1));
            Tau_passk.elb.l = f_passiveTATorques(stiffnessArm, dampingArm, ...
                Xk_nsc(jointi.elb.l*2-1,1), Xk_nsc(jointi.elb.l*2,1));
            Tau_passk.elb.r = f_passiveTATorques(stiffnessArm, dampingArm, ...
                Xk_nsc(jointi.elb.r*2-1,1), Xk_nsc(jointi.elb.r*2,1));        
            Tau_passk.arm = [Tau_passk.sh_flex.l; Tau_passk.sh_add.l; ...
                Tau_passk.sh_rot.l; Tau_passk.sh_flex.r; Tau_passk.sh_add.r; ...
                Tau_passk.sh_rot.r; Tau_passk.elb.l; Tau_passk.elb.r];       

            Tau_passk.mtp.l = f_passiveTATorques(stiffnessMtp, dampingMtp, ...
                Xk_nsc(jointi.mtp.l*2-1,1), Xk_nsc(jointi.mtp.l*2,1));
            Tau_passk.mtp.r = f_passiveTATorques(stiffnessMtp, dampingMtp, ...
                Xk_nsc(jointi.mtp.r*2-1,1), Xk_nsc(jointi.mtp.r*2,1));       
            Tau_passk.mtp.all = [Tau_passk.mtp.l; Tau_passk.mtp.r];

            % Loop over collocation points
            Xk_nsc_end          = D(1)*Xk_nsc;
            FTtildek_nsc_end    = D(1)*FTtildek_nsc;
            ak_end              = D(1)*ak;
            a_ak_end            = D(1)*a_ak;
            a_mtpk_end          = D(1)*a_mtpk;
            for j=1:d
                % Expression for the state derivatives at the collocation point
                xp_nsc          = C(1,j+1)*Xk_nsc;
                FTtildep_nsc    = C(1,j+1)*FTtildek_nsc;
                ap              = C(1,j+1)*ak;
                a_ap            = C(1,j+1)*a_ak;
                a_mtpp          = C(1,j+1)*a_mtpk;
                for r=1:d
                    xp_nsc       = xp_nsc + C(r+1,j+1)*Xkj_nsc{r};
                    FTtildep_nsc = FTtildep_nsc + C(r+1,j+1)*FTtildekj_nsc{r};
                    ap           = ap + C(r+1,j+1)*akj{r};
                    a_ap         = a_ap + C(r+1,j+1)*a_akj{r};
                    a_mtpp       = a_mtpp + C(r+1,j+1)*a_mtpkj{r};
                end 
                % Append collocation equations
                % Dynamic constraints are scaled using the same scale
                % factors as was used to scale the states
                % Activation dynamics (implicit formulation)  
                g       = {g{:}, (h*vAk.*scaling(p).p.vA - ap)./scaling(p).p.a};
                lbg     = [lbg; zeros(NMuscle,1)];
                ubg     = [ubg; zeros(NMuscle,1)]; 
                % Contraction dynamics (implicit formulation)               
                g       = {g{:}, (h*dFTtildek.*scaling(p).p.dFTtilde - ...
                    FTtildep_nsc)./(scaling(p).p.FTtilde')};
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
                    Xkj_nsc{j}(58); Ak_nsc(29); Xkj_nsc{j}(60); Ak_nsc(30); ...
                    Xkj_nsc{j}(62); Ak_nsc(31)];
                g       = {g{:}, (h*xj_nsc - xp_nsc)./(scaling(p).p.QsQdots')};
                lbg     = [lbg; zeros(2*nq.all,1)];
                ubg     = [ubg; zeros(2*nq.all,1)];   
                % Arm activation dynamics (explicit formulation)   
                dadt    = f_ArmActivationDynamics(e_ak,a_akj{j});
                g       = {g{:}, (h*dadt - a_ap)./scaling(p).p.a_a};
                lbg     = [lbg; zeros(nq.arms,1)];
                ubg     = [ubg; zeros(nq.arms,1)]; 
                % Mtp activation dynamics (explicit formulation)
                da_mtpdt    = f_MtpActivationDynamics(e_mtpk,a_mtpkj{j});
                g           = {g{:}, (h*da_mtpdt - a_mtpp)./scaling(p).p.a_mtp};
                lbg         = [lbg; zeros(nq.mtp,1)];
                ubg         = [ubg; zeros(nq.mtp,1)]; 

                % Add contribution to the end state
                Xk_nsc_end = Xk_nsc_end + D(j+1)*Xkj_nsc{j};
                FTtildek_nsc_end = FTtildek_nsc_end + D(j+1)*FTtildekj_nsc{j};
                ak_end = ak_end + D(j+1)*akj{j};  
                a_ak_end = a_ak_end + D(j+1)*a_akj{j};   
                a_mtpk_end = a_mtpk_end + D(j+1)*a_mtpkj{j}; 
                % Add contribution to quadrature function
                J = J + ...
                    W.Qs*B(j+1)*(f_J30(Xk(Qsi(residual_bptyi),1)-...
                        Qs(p).p.allinterpfilt(k+1,residual_bptyi+1)'... 
                        ./scaling(p).p.Qs(residual_bptyi)'))*h +...
                    W.GRF*B(j+1)*(f_J6((Tk(GRFi.all,1)./scaling(p).p.GRF')-...
                        GRF(p).p.val.allinterp(k+1,2:end)'./scaling(p).p.GRF'))*h +...
                    W.GRM*B(j+1)*(f_J6((Tk(GRMi.all,1)./scaling(p).p.GRM')-...
                        GRF(p).p.MorGF.allinterp(k+1,2:end)'./scaling(p).p.GRM'))*h +...
                    W.ID_act*B(j+1)*(f_J23((Tk(residuals_act_bmtpi,1)./scaling(p).p.T(1)')-...
                        ID(p).p.allinterp(k+1,residuals_act_bmtpi+1)'./scaling(p).p.T(1)))*h +...
                    W.a*B(j+1)*(f_J92(akj{j}))*h + ...
                    W.a*B(j+1)*(f_J2(e_mtpk))*h + ...
                    W.a*B(j+1)*(f_J8(e_ak))*h + ...
                    W.u*B(j+1)*(f_J31(Ak))*h +...
                    W.u*B(j+1)*(f_J92(vAk))*h +...
                    W.u*B(j+1)*(f_J92(dFTtildek))*h;
            end                              
            % Add path constraints
            % Pelvis residuals (same as from inverse dynamics)
            g   = {g{:}, Tk(ground_pelvisi,1)};
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
            Ft_trunk_ext    = FTk([mai(8).mus.l,mai(8).mus.r]',1);
            T_trunk_ext     = f_T6(MA.trunk_ext,Ft_trunk_ext);
            g               = {g{:},Tk(jointi.trunk.ext,1) - ...
                (T_trunk_ext + Tau_passk.trunk.ext)};
            lbg             = [lbg; 0];
            ubg             = [ubg; 0];
            % Lumbar bending
            Ft_trunk_ben    = FTk([mai(9).mus.l,mai(9).mus.r]',1);
            T_trunk_ben     = f_T6(MA.trunk_ben,Ft_trunk_ben);
            g               = {g{:},Tk(jointi.trunk.ben,1) - ...
                (T_trunk_ben + Tau_passk.trunk.ben)};
            lbg             = [lbg; 0];
            ubg             = [ubg; 0];
            % Lumbar rotation
            Ft_trunk_rot    = FTk([mai(10).mus.l,mai(10).mus.r]',1);
            T_trunk_rot     = f_T6(MA.trunk_rot,Ft_trunk_rot);
            g               = {g{:},Tk(jointi.trunk.rot,1) - ...
                (T_trunk_rot + Tau_passk.trunk.rot)};
            lbg             = [lbg; 0];
            ubg             = [ubg; 0];
            % Torque-driven joint torques for the arms
            % Arms
            g               = {g{:},Tk(armsi,1)/scaling(p).p.ArmTau - ...
                (a_ak + Tau_passk.arm/scaling(p).p.ArmTau)};
            lbg             = [lbg; zeros(nq.arms,1)];
            ubg             = [ubg; zeros(nq.arms,1)];        
            % Torque-driven joint torques for the mtps
            g               = {g{:},Tk(mtpi,1)/scaling(p).p.MtpTau - ...
                (a_mtpk + Tau_passk.mtp.all/scaling(p).p.MtpTau)};
            lbg             = [lbg; zeros(nq.mtp,1)];
            ubg             = [ubg; zeros(nq.mtp,1)];            
            % Activation dynamics (implicit formulation)
            tact = 0.015;
            tdeact = 0.06;
            act1 = vAk*scaling(p).p.vA + ak./(ones(size(ak,1),1)*tdeact);
            act2 = vAk*scaling(p).p.vA + ak./(ones(size(ak,1),1)*tact);
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
                lbw             = [lbw; bounds(p).p.a.lower'];
                ubw             = [ubw; bounds(p).p.a.upper'];
                w0              = [w0;  guess(p).p.a(k+2,:)'];
                % Muscle-tendon forces
                FTtildek        = MX.sym(['FTtilde_' num2str(k+1)], NMuscle);
                w               = {w{:}, FTtildek};
                lbw             = [lbw; bounds(p).p.FTtilde.lower'];
                ubw             = [ubw; bounds(p).p.FTtilde.upper'];
                w0              = [w0;  guess(p).p.FTtilde(k+2,:)'];    
                % Qs and Qdots
                Xk              = MX.sym(['X_' num2str(k+1)],2*nq.all);
                w               = {w{:}, Xk};
                lbw             = [lbw; bounds(p).p.QsQdots.lower'];
                ubw             = [ubw; bounds(p).p.QsQdots.upper']; 
                w0              = [w0;  guess(p).p.QsQdots(k+2,:)'];
                % Arm activations
                a_ak            = MX.sym(['a_a_' num2str(k+1)], nq.arms);
                w               = {w{:}, a_ak};
                lbw             = [lbw; bounds(p).p.a_a.lower'];
                ubw             = [ubw; bounds(p).p.a_a.upper'];
                w0              = [w0;  guess(p).p.a_a(k+2,:)'];
                % Mtp activations
                a_mtpk          = MX.sym(['a_mtp_' num2str(k+1)], nq.mtp);
                w               = {w{:}, a_mtpk};
                lbw             = [lbw; bounds(p).p.a_mtp.lower'];
                ubw             = [ubw; bounds(p).p.a_mtp.upper'];
                w0              = [w0;  guess(p).p.a_mtp(k+2,:)'];
            else
                % Muscle activations
                ak              = MX.sym(['a_' num2str(k+1)], NMuscle);
                w               = {w{:}, ak};
                lbw             = [lbw; bounds(p).p.a.lower'];
                ubw             = [ubw; bounds(p).p.a.upper'];
                w0              = [w0;  guess(p).p.a(end,:)'];
                % Muscle-tendon forces
                FTtildek        = MX.sym(['FTtilde_' num2str(k+1)], NMuscle);
                w               = {w{:}, FTtildek};
                lbw             = [lbw; bounds(p).p.FTtilde.lower'];
                ubw             = [ubw; bounds(p).p.FTtilde.upper'];
                w0              = [w0;  guess(p).p.FTtilde(end,:)'];    
                % Qs and Qdots
                Xk              = MX.sym(['X_' num2str(k+1)],2*nq.all);
                w               = {w{:}, Xk};
                lbw             = [lbw; bounds(p).p.QsQdots.lower'];
                ubw             = [ubw; bounds(p).p.QsQdots.upper']; 
                w0              = [w0;  guess(p).p.QsQdots(end,:)'];
                % Arm activations
                a_ak            = MX.sym(['a_a_' num2str(k+1)], nq.arms);
                w               = {w{:}, a_ak};
                lbw             = [lbw; bounds(p).p.a_a.lower'];
                ubw             = [ubw; bounds(p).p.a_a.upper'];
                w0              = [w0;  guess(p).p.a_a(end,:)'];
                % Mtp activations
                a_mtpk          = MX.sym(['a_mtp_' num2str(k+1)], nq.mtp);
                w               = {w{:}, a_mtpk};
                lbw             = [lbw; bounds(p).p.a_mtp.lower'];
                ubw             = [ubw; bounds(p).p.a_mtp.upper'];
                w0              = [w0;  guess(p).p.a_mtp(end,:)'];
            end
            % Rescale variables to impose equality constraints
            Xk_end = (Xk_nsc_end)./scaling(p).p.QsQdots';
            FTtildek_end = (FTtildek_nsc_end)./scaling(p).p.FTtilde';
            % Add equality constraints (next interval starts with end values of 
            % states from previous interval)
            g   = {g{:}, Xk_end-Xk, FTtildek_end-FTtildek, ...
                ak_end-ak, a_ak_end-a_ak, a_mtpk_end-a_mtpk};
            lbg = [lbg; zeros(2*nq.all + NMuscle + NMuscle + nq.arms + nq.mtp,1)];
            ubg = [ubg; zeros(2*nq.all + NMuscle + NMuscle + nq.arms + nq.mtp,1)];    

        end       
        % Periodicity pelvis_ty (Qs and Qdots)
        pelvis_tyi = 2*jointi.pelvis.ty-1:2*jointi.pelvis.ty;
        g   = {g{:}, Xk_end(pelvis_tyi)-X0(pelvis_tyi,1)};
        lbg = [lbg; zeros(length(pelvis_tyi),1)];
        ubg = [ubg; zeros(length(pelvis_tyi),1)];    
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
    stats = solver.stats();
    % Save results and setup
    save([pathresults,'/',namescript,'/w',savename],'w_opt');
    save([pathresults,'/',namescript,'/stats',savename],'stats');
end

%% Analyze results
if analyseResults
    %% Load results
    if loadResults
        p = mfilename('fullpath');
        [~,namescript,~] = fileparts(p);
        pathresults = [pathRepo,'/Results'];
        load([pathresults,'/',namescript,'/w',savename]);
        load([pathresults,'/',namescript,'/stats',savename]);
    end
    
    %% Extract results
    % All optimized design variables are saved in a single column vector      
    % Number of design variables    
    NControls = NMuscle+NMuscle+nq.all+nq.arms+nq.mtp;
    NStates = NMuscle+NMuscle+2*nq.all+nq.arms+nq.mtp;
    NParameters = np;
    % In the loop
    Nwl = NControls+d*NStates+NStates;
    Nw = NParameters+NStates+N*Nwl;
    % Before the variable corresponding to the first collocation point
    Nwm = NParameters+NStates+NControls;
    % Here we extract the results and re-organize them for analysis  
    % Static parameters
    paramsCM_opt = w_opt(1:NParameters);  
    N4p = 0;
    for p = 1:length(trials)     
        % Mesh points
        % Muscle activations and muscle-tendon forces
        a_opt(p).p = zeros(N+1,NMuscle);
        FTtilde_opt(p).p = zeros(N+1,NMuscle);
        for i = 1:NMuscle
            a_opt(p).p(:,i) = w_opt(NParameters+N4p+i:Nwl:Nw);
            FTtilde_opt(p).p(:,i) = w_opt(NParameters+N4p+NMuscle+i:Nwl:Nw);
        end
        % Qs and Qdots
        q_opt(p).p = zeros(N+1,nq.all);
        qdot_opt(p).p = zeros(N+1,nq.all);
        count = 0;
        for i = 1:2:2*nq.all
            count = count +1;
            q_opt(p).p(:,count) = w_opt(NParameters+N4p+NMuscle+NMuscle+i:Nwl:Nw);
            qdot_opt(p).p(:,count) = w_opt(NParameters+N4p+NMuscle+NMuscle+i+1:Nwl:Nw);
        end
        % Arm activations
        a_a_opt(p).p = zeros(N+1,nq.arms);
        for i = 1:nq.arms
            a_a_opt(p).p(:,i) = w_opt(NParameters+N4p+NMuscle+NMuscle+2*nq.all+i:Nwl:Nw);
        end
        % Mtp activations
        a_mtp_opt(p).p = zeros(N+1,nq.mtp);
        for i = 1:nq.mtp
            a_mtp_opt(p).p(:,i) = w_opt(NParameters+N4p+NMuscle+NMuscle+2*nq.all+...
                nq.arms+i:Nwl:Nw);
        end
        % Time derivative of muscle activations and muscle-tendon forces
        vA_opt(p).p = zeros(N,NMuscle);
        dFTtilde_opt(p).p = zeros(N,NMuscle);
        for i = 1:NMuscle
            vA_opt(p).p(:,i) = w_opt(NParameters+N4p+NStates+i:Nwl:Nw);
            dFTtilde_opt(p).p(:,i) = w_opt(NParameters+N4p+NStates+i:Nwl:Nw);
        end
        % Time derivative of joint velocities
        qdotdot_opt(p).p = zeros(N,nq.all);
        for i = 1:nq.all
            qdotdot_opt(p).p(:,i) = w_opt(NParameters+N4p+NStates+2*NMuscle+i:Nwl:Nw);
        end
        % Arm excitations
        e_a_opt(p).p = zeros(N,nq.arms);
        for i = 1:nq.arms
            e_a_opt(p).p(:,i) = w_opt(NParameters+N4p+NStates+2*NMuscle+nq.all+i:Nwl:Nw); 
        end 
        % Mtp excitations
        e_mtp_opt(p).p = zeros(N,nq.mtp);
        for i = 1:nq.arms
            e_mtp_opt(p).p(:,i) = w_opt(NParameters+N4p+NStates+2*NMuscle+nq.all+...
                nq.arms+i:Nwl:Nw); 
        end    
        N4p = Nw - NParameters;
        Nw = Nw+NStates+N*Nwl;
%         % Collocation points
%         % Muscle activations
%         a_opt_ext(p).p=zeros(N*(d+1)+1,NMuscle);
%         a_opt_ext(p).p(1:(d+1):end,:)= a_opt(p).p;
%         for nmusi=1:NMuscle
%             a_opt_ext(p).p(2:(d+1):end,nmusi) = w_opt(Nwm+nmusi:Nwl:Nw);
%             a_opt_ext(p).p(3:(d+1):end,nmusi) = ...
%                 w_opt(Nwm+NMuscle+nmusi:Nwl:Nw);
%             a_opt_ext(p).p(4:(d+1):end,nmusi) = ...
%                 w_opt(Nwm+NMuscle+NMuscle+nmusi:Nwl:Nw);
%         end
%         % Muscle-tendon forces
%         FTtilde_opt_ext(p).p=zeros(N*(d+1)+1,NMuscle);
%         FTtilde_opt_ext(p).p(1:(d+1):end,:)= FTtilde_opt(p).p;
%         for nmusi=1:NMuscle
%             FTtilde_opt_ext(p).p(2:(d+1):end,nmusi) = ...
%                 w_opt(Nwm+d*NMuscle+nmusi:Nwl:Nw);
%             FTtilde_opt_ext(p).p(3:(d+1):end,nmusi) = ...
%                 w_opt(Nwm+d*NMuscle+NMuscle+nmusi:Nwl:Nw);
%             FTtilde_opt_ext(p).p(4:(d+1):end,nmusi) = ...
%                 w_opt(Nwm+d*NMuscle+NMuscle+NMuscle+nmusi:Nwl:Nw);
%         end
%         % Qs and Qdots
%         q_opt_ext(p).p=zeros(N*(d+1)+1,nq.all);
%         q_opt_ext(p).p(1:(d+1):end,:)= q_opt(p).p;
%         q_dot_opt_ext=zeros(N*(d+1)+1,nq.all);
%         q_dot_opt_ext(1:(d+1):end,:)= qdot_opt(p).p;
%         nqi_col = 1:2:2*nq.all;
%         for nqi=1:nq.all
%             nqi_q = nqi_col(nqi);
%             q_opt_ext(p).p(2:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
%                 d*NMuscle+nqi_q:Nwl:Nw);   
%             q_opt_ext(p).p(3:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
%                 d*NMuscle+2*nq.all+nqi_q:Nwl:Nw);  
%             q_opt_ext(p).p(4:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
%                 d*NMuscle+2*nq.all+2*nq.all+nqi_q:Nwl:Nw);  
%             q_dot_opt_ext(2:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
%                 d*NMuscle+nqi_q+1:Nwl:Nw);   
%             q_dot_opt_ext(3:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
%                 d*NMuscle+2*nq.all+nqi_q+1:Nwl:Nw);  
%             q_dot_opt_ext(4:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
%                 d*NMuscle+2*nq.all+2*nq.all+nqi_q+1:Nwl:Nw);
%         end
%         % Arm activations
%         a_a_opt_ext(p).p=zeros(N*(d+1)+1,nq.arms);
%         a_a_opt_ext(p).p(1:(d+1):end,:)= a_a_opt(p).p;
%         for nmusi=1:nq.arms
%             a_a_opt_ext(p).p(2:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle+...
%                 d*NMuscle+d*2*nq.all+nmusi:Nwl:Nw);
%             a_a_opt_ext(p).p(3:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle+...
%                 d*NMuscle+d*2*nq.all+nq.arms+nmusi:Nwl:Nw);
%             a_a_opt_ext(p).p(4:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle+...
%                 d*NMuscle+d*2*nq.all+nq.arms+nq.arms+nmusi:Nwl:Nw);
%         end
%         % Mtp activations
%         a_mtp_opt_ext(p).p=zeros(N*(d+1)+1,nq.mtp);
%         a_mtp_opt_ext(p).p(1:(d+1):end,:)= a_mtp_opt(p).p;
%         for nmusi=1:nq.mtp
%             a_mtp_opt_ext(p).p(2:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle+...
%                 d*NMuscle+d*2*nq.all+nq.arms+nmusi:Nwl:Nw);
%             a_mtp_opt_ext(p).p(3:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle+...
%                 d*NMuscle+d*2*nq.all+nq.arms+nq.mtp+nmusi:Nwl:Nw);
%             a_mtp_opt_ext(p).p(4:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle+...
%                 d*NMuscle+d*2*nq.all+nq.arms+nq.mtp+nq.mtp+nmusi:Nwl:Nw);
%         end
    end
    
    %% Unscale results
    % Parameters
    paramsCM_opt_nsc = f_nscnp(paramsCM_opt,scaling(1).p.params.v,scaling(1).p.params.r);
    for p = 1:length(trials)   
        % States at mesh points
        % Qs (1:N-1)
        q_opt_unsc(p).p.rad = q_opt(p).p(1:end-1,:).*repmat(...
        scaling(p).p.Qs,size(q_opt(p).p(1:end-1,:),1),1); 
        % Convert in degrees
        q_opt_unsc(p).p.deg = q_opt_unsc(p).p.rad;
        q_opt_unsc(p).p.deg(:,[1:3,7:end]) = q_opt_unsc(p).p.deg(:,[1:3,7:end]).*180/pi;
        % Qs (1:N)
        q_opt_unsc_all(p).p.rad = q_opt(p).p(1:end,:).*repmat(...
            scaling(p).p.Qs,size(q_opt(p).p(1:end,:),1),1); 
        % Convert in degrees
        q_opt_unsc_all(p).p.deg = q_opt_unsc_all(p).p.rad;
        q_opt_unsc_all(p).p.deg(:,[1:3,7:end]) = ...
            q_opt_unsc_all(p).p.deg(:,[1:3,7:end]).*180/pi;
        % Qdots (1:N-1)
        qdot_opt_unsc(p).p.rad = qdot_opt(p).p(1:end-1,:).*repmat(...
            scaling(p).p.Qdots,size(qdot_opt(p).p(1:end-1,:),1),1);
        % Convert in degrees
        qdot_opt_unsc(p).p.deg = qdot_opt_unsc(p).p.rad;
        qdot_opt_unsc(p).p.deg(:,[1:3,7:end]) = ...
            qdot_opt_unsc(p).p.deg(:,[1:3,7:end]).*180/pi;
        % Qdots (1:N)
        qdot_opt_unsc_all(p).p.rad = qdot_opt(p).p(1:end,:).*repmat(...
            scaling(p).p.Qdots,size(qdot_opt(p).p(1:end,:),1),1); 
        % Muscle activations
        a_opt_unsc(p).p = a_opt(p).p(1:end-1,:).*repmat(...
            scaling(p).p.a,size(a_opt(p).p(1:end-1,:),1),size(a_opt(p).p,2));
        % Muscle-tendon forces
        FTtilde_opt_unsc(p).p = FTtilde_opt(p).p(1:end-1,:).*repmat(...
            scaling(p).p.FTtilde,size(FTtilde_opt(p).p(1:end-1,:),1),1);
        % Arm activations
        a_a_opt_unsc(p).p = a_a_opt(p).p(1:end-1,:);
        % Mtp activations
        a_mtp_opt_unsc(p).p = a_mtp_opt(p).p(1:end-1,:);
        % Controls at mesh points
        % Time derivative of Qdots
        qdotdot_opt_unsc(p).p.rad = ...
            qdotdot_opt(p).p.*repmat(scaling(p).p.Qdotdots,size(qdotdot_opt(p).p,1),1);
        % Convert in degrees
        qdotdot_opt_unsc(p).p.deg = qdotdot_opt_unsc(p).p.rad;
        qdotdot_opt_unsc(p).p.deg(:,[1:3,7:end]) = ...
            qdotdot_opt_unsc(p).p.deg(:,[1:3,7:end]).*180/pi;
        % Time derivative of muscle activations (states)
        vA_opt_unsc(p).p = vA_opt(p).p.*repmat(scaling(p).p.vA,size(vA_opt(p).p,1),size(vA_opt(p).p,2));
        tact = 0.015;
        tdeact = 0.06;
        % Get muscle excitations from time derivative of muscle activations
        e_opt_unsc(p).p = computeExcitationRaasch(a_opt_unsc(p).p,vA_opt_unsc(p).p,...
            ones(1,NMuscle)*tdeact,ones(1,NMuscle)*tact);
        % Time derivative of muscle-tendon forces
        dFTtilde_opt_unsc(p).p = dFTtilde_opt(p).p.*repmat(...
            scaling(p).p.dFTtilde,size(dFTtilde_opt(p).p,1),size(dFTtilde_opt(p).p,2));
        % Arm excitations
        e_a_opt_unsc(p).p = e_a_opt(p).p.*repmat(scaling(p).p.e_a,size(e_a_opt(p).p,1),...
            size(e_a_opt(p).p,2));
        % Mtp excitations
        e_mtp_opt_unsc(p).p = e_mtp_opt(p).p.*repmat(scaling(p).p.e_mtp,size(e_mtp_opt(p).p,1),...
            size(e_mtp_opt(p).p,2));
    end
    
    %% Time grid    
    % Mesh points
    for p = 1:length(trials) 
        tgrid(p).p = linspace(time_opt(p).p(1),time_opt(p).p(end),N+1);
        dtime = zeros(1,d+1);
        for i=1:4
            dtime(i)=tau_root(i)*((time_opt(p).p(end)-time_opt(p).p(1))/N);
        end
        % Mesh points and collocation points
        tgrid_ext(p).p = zeros(1,(d+1)*N+1);
        for i=1:N
            tgrid_ext(p).p(((i-1)*4+1):1:i*4)=tgrid(p).p(i)+dtime;
        end
        tgrid_ext(p).p(end)=time_opt(p).p(end); 
    end
 
    %% Joint torques, ground reaction forces and moments at opt solution
    for p = 1:length(trials)    
        Xk_Qs_Qdots_opt             = zeros(N,2*size(q_opt_unsc(p).p.rad,2));
        Xk_Qs_Qdots_opt(:,1:2:end)  = q_opt_unsc(p).p.rad;
        Xk_Qs_Qdots_opt(:,2:2:end)  = qdot_opt_unsc(p).p.rad;
        Xk_Qdotdots_opt             = qdotdot_opt_unsc(p).p.rad;  
        out1_res_opt                = zeros(N,toes.TR.T.r(end));
        out2_res_opt                = zeros(N,GRMi.l(end));
        Tau_pass_opt_mtp            = zeros(N,nq.mtp);
        Tau_pass_opt_arm            = zeros(N,nq.arms);
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
            if cs == 5 || cs == 6
                locSphere_s5_l_opt = [paramsCM_opt_nsc(loci.s5.r.x),...
                    loc.s5.y,-paramsCM_opt_nsc(loci.s5.r.z)]';
                locSphere_s5_r_opt = [paramsCM_opt_nsc(loci.s5.r.x),...
                    loc.s5.y,paramsCM_opt_nsc(loci.s5.r.z)]';
                if cs == 6
                    locSphere_s6_l_opt = [paramsCM_opt_nsc(loci.s6.r.x),...
                        loc.s6.y,-paramsCM_opt_nsc(loci.s6.r.z)]';
                    locSphere_s6_r_opt = [paramsCM_opt_nsc(loci.s6.r.x),...
                        loc.s6.y,paramsCM_opt_nsc(loci.s6.r.z)]';   
                end
            end
            % Compute contact forces
            force_s1_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s1),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s1_l_opt,...
                out1_res_opt(i,calcn.l.pos),out1_res_opt(i,calcn.l.v_lin),...
                out1_res_opt(i,calcn.l.omega),out1_res_opt(i,calcn.TR.R.l),...
                out1_res_opt(i,calcn.TR.T.l));
            force_s1_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s1),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s1_r_opt,...
                out1_res_opt(i,calcn.r.pos),out1_res_opt(i,calcn.r.v_lin),...
                out1_res_opt(i,calcn.r.omega),out1_res_opt(i,calcn.TR.R.r),...
                out1_res_opt(i,calcn.TR.T.r));
            force_s2_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s2),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s2_l_opt,...
                out1_res_opt(i,calcn.l.pos),out1_res_opt(i,calcn.l.v_lin),...
                out1_res_opt(i,calcn.l.omega),out1_res_opt(i,calcn.TR.R.l),...
                out1_res_opt(i,calcn.TR.T.l)); 
            force_s2_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s2),...
                dissipation,normal,transitionVelo,staticFriction,...
                dynamicFriction,viscousFriction,locSphere_s2_r_opt,...
                out1_res_opt(i,calcn.r.pos),out1_res_opt(i,calcn.r.v_lin),...
                out1_res_opt(i,calcn.r.omega),out1_res_opt(i,calcn.TR.R.r),...
                out1_res_opt(i,calcn.TR.T.r));  
            if (cs == 4 && csc == 2) || (cs == 4 && csc == 3) % Sphere 3 => toes  
                force_s3_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s3),...
                    dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,locSphere_s3_l_opt,...
                    out1_res_opt(i,toes.l.pos),out1_res_opt(i,toes.l.v_lin),...
                    out1_res_opt(i,toes.l.omega),out1_res_opt(i,toes.TR.R.l),...
                    out1_res_opt(i,toes.TR.T.l));  
                force_s3_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s3),...
                    dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,locSphere_s3_r_opt,...
                    out1_res_opt(i,toes.r.pos),out1_res_opt(i,toes.r.v_lin),...
                    out1_res_opt(i,toes.r.omega),out1_res_opt(i,toes.TR.R.r),...
                    out1_res_opt(i,toes.TR.T.r));
            else % Sphere 3 => calcn
                force_s3_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s3),...
                    dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,locSphere_s3_l_opt,...
                    out1_res_opt(i,calcn.l.pos),out1_res_opt(i,calcn.l.v_lin),...
                    out1_res_opt(i,calcn.l.omega),out1_res_opt(i,calcn.TR.R.l),...
                    out1_res_opt(i,calcn.TR.T.l));
                force_s3_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s3),...
                    dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,locSphere_s3_r_opt,...
                    out1_res_opt(i,calcn.r.pos),out1_res_opt(i,calcn.r.v_lin),...
                    out1_res_opt(i,calcn.r.omega),out1_res_opt(i,calcn.TR.R.r),...
                    out1_res_opt(i,calcn.TR.T.r));            
            end
            if (cs == 6 && csc == 4) || (cs == 6 && csc == 5) || (cs == 6 && csc == 6) % Sphere 4 => calcn
                force_s4_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s4),...
                    dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,locSphere_s4_l_opt,...
                    out1_res_opt(i,calcn.l.pos),out1_res_opt(i,calcn.l.v_lin),...
                    out1_res_opt(i,calcn.l.omega),out1_res_opt(i,calcn.TR.R.l),...
                    out1_res_opt(i,calcn.TR.T.l));  
                force_s4_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s4),...
                    dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,locSphere_s4_r_opt,...
                    out1_res_opt(i,calcn.r.pos),out1_res_opt(i,calcn.r.v_lin),...
                    out1_res_opt(i,calcn.r.omega),out1_res_opt(i,calcn.TR.R.r),...
                    out1_res_opt(i,calcn.TR.T.r));                 
            else % Sphere 4 => toes
                force_s4_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s4),...
                    dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,locSphere_s4_l_opt,...
                    out1_res_opt(i,toes.l.pos),out1_res_opt(i,toes.l.v_lin),...
                    out1_res_opt(i,toes.l.omega),out1_res_opt(i,toes.TR.R.l),...
                    out1_res_opt(i,toes.TR.T.l));  
                force_s4_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s4),...
                    dissipation,normal,transitionVelo,staticFriction,...
                    dynamicFriction,viscousFriction,locSphere_s4_r_opt,...
                    out1_res_opt(i,toes.r.pos),out1_res_opt(i,toes.r.v_lin),...
                    out1_res_opt(i,toes.r.omega),out1_res_opt(i,toes.TR.R.r),...
                    out1_res_opt(i,toes.TR.T.r));                  
            end           
            if cs == 5 || cs == 6
                if csc == 1 % Sphere 5 => calcn
                    force_s5_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s5),...
                        dissipation,normal,transitionVelo,staticFriction,...
                        dynamicFriction,viscousFriction,locSphere_s5_r_opt,...
                        out1_res_opt(i,calcn.r.pos),out1_res_opt(i,calcn.r.v_lin),...
                        out1_res_opt(i,calcn.r.omega),out1_res_opt(i,calcn.TR.R.r),...
                        out1_res_opt(i,calcn.TR.T.r));
                    force_s5_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s5),...
                        dissipation,normal,transitionVelo,staticFriction,...
                        dynamicFriction,viscousFriction,locSphere_s5_l_opt,...
                        out1_res_opt(i,calcn.l.pos),out1_res_opt(i,calcn.l.v_lin),...
                        out1_res_opt(i,calcn.l.omega),out1_res_opt(i,calcn.TR.R.l),...
                        out1_res_opt(i,calcn.TR.T.l));
                elseif csc == 2 || csc == 3 || csc == 4 || csc == 5 || csc == 6 % Sphere 5 => toes
                    force_s5_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s5),...
                        dissipation,normal,transitionVelo,staticFriction,...
                        dynamicFriction,viscousFriction,locSphere_s5_r_opt,...
                        out1_res_opt(i,toes.r.pos),out1_res_opt(i,toes.r.v_lin),...
                        out1_res_opt(i,toes.r.omega),out1_res_opt(i,toes.TR.R.r),...
                        out1_res_opt(i,toes.TR.T.r));
                    force_s5_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s5),...
                        dissipation,normal,transitionVelo,staticFriction,...
                        dynamicFriction,viscousFriction,locSphere_s5_l_opt,...
                        out1_res_opt(i,toes.l.pos),out1_res_opt(i,toes.l.v_lin),...
                        out1_res_opt(i,toes.l.omega),out1_res_opt(i,toes.TR.R.l),...
                        out1_res_opt(i,toes.TR.T.l));                    
                end
                if cs == 6
                    force_s6_r_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s6),...
                        dissipation,normal,transitionVelo,staticFriction,...
                        dynamicFriction,viscousFriction,locSphere_s6_r_opt,...
                        out1_res_opt(i,toes.r.pos),out1_res_opt(i,toes.r.v_lin),...
                        out1_res_opt(i,toes.r.omega),out1_res_opt(i,toes.TR.R.r),...
                        out1_res_opt(i,toes.TR.T.r));  
                    force_s6_l_opt = f_contactForce(stif,paramsCM_opt_nsc(radi.s6),...
                        dissipation,normal,transitionVelo,staticFriction,...
                        dynamicFriction,viscousFriction,locSphere_s6_l_opt,...
                        out1_res_opt(i,toes.l.pos),out1_res_opt(i,toes.l.v_lin),...
                        out1_res_opt(i,toes.l.omega),out1_res_opt(i,toes.TR.R.l),...
                        out1_res_opt(i,toes.TR.T.l));    
                end
            end
            if cs == 4
                in_F2_opt = [force_s1_l_opt,force_s2_l_opt,force_s3_l_opt,...
                    force_s4_l_opt,force_s1_r_opt,...
                    force_s2_r_opt,force_s3_r_opt,force_s4_r_opt,...
                    paramsCM_opt_nsc(loci.s1.r.x),paramsCM_opt_nsc(loci.s1.r.z),...
                    paramsCM_opt_nsc(loci.s2.r.x),paramsCM_opt_nsc(loci.s2.r.z),...
                    paramsCM_opt_nsc(loci.s3.r.x),paramsCM_opt_nsc(loci.s3.r.z),...
                    paramsCM_opt_nsc(loci.s4.r.x),paramsCM_opt_nsc(loci.s4.r.z),...
                    paramsCM_opt_nsc(radi.s1:radi.s4)'];   
            elseif cs == 5
                in_F2_opt = [force_s1_l_opt,force_s2_l_opt,force_s3_l_opt,...
                    force_s4_l_opt,force_s5_l_opt,force_s1_r_opt,...
                    force_s2_r_opt,force_s3_r_opt,force_s4_r_opt,force_s5_r_opt,...
                    paramsCM_opt_nsc(loci.s1.r.x),paramsCM_opt_nsc(loci.s1.r.z),...
                    paramsCM_opt_nsc(loci.s2.r.x),paramsCM_opt_nsc(loci.s2.r.z),...
                    paramsCM_opt_nsc(loci.s3.r.x),paramsCM_opt_nsc(loci.s3.r.z),...
                    paramsCM_opt_nsc(loci.s4.r.x),paramsCM_opt_nsc(loci.s4.r.z),...
                    paramsCM_opt_nsc(loci.s5.r.x),paramsCM_opt_nsc(loci.s5.r.z),...
                    paramsCM_opt_nsc(radi.s1:radi.s5)'];   
            elseif cs == 6
                in_F2_opt = [force_s1_l_opt,force_s2_l_opt,force_s3_l_opt,...
                    force_s4_l_opt,force_s5_l_opt,force_s6_l_opt,force_s1_r_opt,...
                    force_s2_r_opt,force_s3_r_opt,force_s4_r_opt,force_s5_r_opt,...
                    force_s6_r_opt,...
                    paramsCM_opt_nsc(loci.s1.r.x),paramsCM_opt_nsc(loci.s1.r.z),...
                    paramsCM_opt_nsc(loci.s2.r.x),paramsCM_opt_nsc(loci.s2.r.z),...
                    paramsCM_opt_nsc(loci.s3.r.x),paramsCM_opt_nsc(loci.s3.r.z),...
                    paramsCM_opt_nsc(loci.s4.r.x),paramsCM_opt_nsc(loci.s4.r.z),...
                    paramsCM_opt_nsc(loci.s5.r.x),paramsCM_opt_nsc(loci.s5.r.z),...
                    paramsCM_opt_nsc(loci.s6.r.x),paramsCM_opt_nsc(loci.s6.r.z),...
                    paramsCM_opt_nsc(radi.s1:radi.s6)'];   
            end
            out2_res = ...
                F2([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)';in_F2_opt']);        
            out2_res_opt(i,:) = full(out2_res); 
            
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

            Tau_pass_opt.sh_flex.l = f_passiveTATorques(stiffnessArm, dampingArm,...
                Xk_Qs_Qdots_opt(i,jointi.sh_flex.l*2-1), ...
                Xk_Qs_Qdots_opt(i,jointi.sh_flex.l*2));
            Tau_pass_opt.sh_add.l = f_passiveTATorques(stiffnessArm, dampingArm,...
                Xk_Qs_Qdots_opt(i,jointi.sh_add.l*2-1), ...
                Xk_Qs_Qdots_opt(i,jointi.sh_add.l*2));
            Tau_pass_opt.sh_rot.l = f_passiveTATorques(stiffnessArm, dampingArm,...
                Xk_Qs_Qdots_opt(i,jointi.sh_rot.l*2-1), ...
                Xk_Qs_Qdots_opt(i,jointi.sh_rot.l*2));
            Tau_pass_opt.sh_flex.r = f_passiveTATorques(stiffnessArm, dampingArm,...
                Xk_Qs_Qdots_opt(i,jointi.sh_flex.r*2-1), ...
                Xk_Qs_Qdots_opt(i,jointi.sh_flex.r*2));
            Tau_pass_opt.sh_add.r = f_passiveTATorques(stiffnessArm, dampingArm,...
                Xk_Qs_Qdots_opt(i,jointi.sh_add.r*2-1), ...
                Xk_Qs_Qdots_opt(i,jointi.sh_add.r*2));
            Tau_pass_opt.sh_rot.r = f_passiveTATorques(stiffnessArm, dampingArm,...
                Xk_Qs_Qdots_opt(i,jointi.sh_rot.r*2-1), ...
                Xk_Qs_Qdots_opt(i,jointi.sh_rot.r*2));
            Tau_pass_opt.elb.l = f_passiveTATorques(stiffnessArm, dampingArm, ...
                Xk_Qs_Qdots_opt(i,jointi.elb.l*2-1), ...
                Xk_Qs_Qdots_opt(i,jointi.elb.l*2));
            Tau_pass_opt.elb.r = f_passiveTATorques(stiffnessArm, dampingArm, ...
                Xk_Qs_Qdots_opt(i,jointi.elb.r*2-1), ...
                Xk_Qs_Qdots_opt(i,jointi.elb.r*2));
            
            Tau_pass_opt_arm(i,:) = full([Tau_pass_opt.sh_flex.l; Tau_pass_opt.sh_add.l; ...
                Tau_pass_opt.sh_rot.l; Tau_pass_opt.sh_flex.r; Tau_pass_opt.sh_add.r; ...
                Tau_pass_opt.sh_rot.r; Tau_pass_opt.elb.l; Tau_pass_opt.elb.r]);    

            Tau_pass_opt.mtp.l = f_passiveTATorques(stiffnessMtp, dampingMtp, ...
                Xk_Qs_Qdots_opt(i,jointi.mtp.l*2-1), ...
                Xk_Qs_Qdots_opt(i,jointi.mtp.l*2));
            Tau_pass_opt.mtp.r = f_passiveTATorques(stiffnessMtp, dampingMtp, ...
                Xk_Qs_Qdots_opt(i,jointi.mtp.r*2-1), ...
                Xk_Qs_Qdots_opt(i,jointi.mtp.r*2));  
            
            Tau_pass_opt_mtp(i,:) = full([Tau_pass_opt.mtp.l,Tau_pass_opt.mtp.r]);

            Tau_pass_opt_all(p).p(i,:) = full([Tau_pass_opt.hip.flex.l,...
                Tau_pass_opt.hip.add.l,Tau_pass_opt.hip.rot.l,...             
                Tau_pass_opt.hip.flex.r,Tau_pass_opt.hip.add.r,...
                Tau_pass_opt.hip.rot.r,Tau_pass_opt.knee.l,...
                Tau_pass_opt.knee.r,Tau_pass_opt.ankle.l,...
                Tau_pass_opt.ankle.r,Tau_pass_opt.subt.l,...
                Tau_pass_opt.subt.r,Tau_pass_opt.mtp.l,...
                Tau_pass_opt.mtp.r,Tau_pass_opt.trunk.ext,...
                Tau_pass_opt.trunk.ben,Tau_pass_opt.trunk.rot,...
                Tau_pass_opt.sh_flex.l,Tau_pass_opt.sh_add.l,...
                Tau_pass_opt.sh_rot.l,Tau_pass_opt.sh_flex.r,...
                Tau_pass_opt.sh_add.r,Tau_pass_opt.sh_rot.r,...
                Tau_pass_opt.elb.l,Tau_pass_opt.elb.r]);
            
        end
        % Optimal joint torques, ground reaction forces and moments
        Tauk_out(p).p        = out2_res_opt(:,residualsi);
        GRF_opt_unsc(p).p    = out2_res_opt(:,GRFi.all);
        GRM_opt_unsc(p).p    = out2_res_opt(:,GRMi.all);
        % assertArmTmax should be 0
        assertArmTmax = max(max(abs(out2_res_opt(:,armsi)/scaling(p).p.ArmTau - ...
            (a_a_opt_unsc(p).p + Tau_pass_opt_arm/scaling(p).p.ArmTau))));
        assertMtpTmax = max(max(abs(out2_res_opt(:,mtpi)/scaling(p).p.MtpTau - ...
            (a_mtp_opt_unsc(p).p + Tau_pass_opt_mtp/scaling(p).p.MtpTau))));
        if assertArmTmax > 10^(-tol_ipopt)
            disp('error in arm torques')
        end
        if assertMtpTmax > 10^(-tol_ipopt)
            disp('error in mtp torques')
        end      
    end
    
    %% Create .mot file for OpenSim GUI
    if writeIKmotion  
        pathOpenSim = [pathRepo,'/OpenSim'];
        addpath(genpath(pathOpenSim));
        for p = 1:length(trials)             
            q_opt_GUI = zeros(N+1,1+nq.all+2);
            q_opt_GUI(:,1) = tgrid(p).p';
            q_opt_GUI(:,2:nq.all+1)  = q_opt_unsc_all(p).p.deg;
            q_opt_GUI(:,end-1:end) = 1.51*180/pi*ones(N+1,2);% pro_sup (locked)
            JointAngle.labels = {'time','pelvis_tilt','pelvis_list',...
                'pelvis_rotation','pelvis_tx','pelvis_ty','pelvis_tz',...
                'hip_flexion_l','hip_adduction_l','hip_rotation_l',...
                'hip_flexion_r','hip_adduction_r','hip_rotation_r',...
                'knee_angle_l','knee_angle_r','ankle_angle_l',...
                'ankle_angle_r','subtalar_angle_l','subtalar_angle_r',...
                'mtp_angle_l','mtp_angle_r', ...
                'lumbar_extension','lumbar_bending','lumbar_rotation',...
                'arm_flex_l','arm_add_l','arm_rot_l',...
                'arm_flex_r','arm_add_r','arm_rot_r',...
                'elbow_flex_l','elbow_flex_r',...
                'pro_sup_l','pro_sup_r'};
            % Combine data joint angles and muscle activations
            JointAngleMuscleAct.data = [q_opt_GUI,a_opt(p).p];
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
                    '/IK',savename,'_',trials{p},'.mot'];
            write_motionFile(JointAngleMuscleAct, filenameJointAngles)
        end
    end   
           
    %% Save results       
    if saveResults
        if (exist([pathresults,'/',namescript,'/Results_tracking.mat'],...
            'file')==2) 
            load([pathresults,'/',namescript,'/Results_tracking.mat']);
        else
            Results_tracking = struct('ww',[]);
        end
    % Structure results
    for p = 1:length(trials) 
        Results_tracking(ww).ww(p).Qs_opt = q_opt_unsc(p).p.deg;
        Results_tracking(ww).ww(p).Acts_opt = a_opt_unsc(p).p;
        Results_tracking(ww).ww(p).Ts_opt = Tauk_out(p).p;
        Results_tracking(ww).ww(p).GRFs_opt = GRF_opt_unsc(p).p;
        Results_tracking(ww).ww(p).GRMs_opt = GRM_opt_unsc(p).p; 
        Results_tracking(ww).ww(p).passT = Tau_pass_opt_all(p).p;      
        Results_tracking(ww).ww(p).aMTP = a_mtp_opt_unsc(p).p;
        Results_tracking(ww).ww(p).aArm = a_a_opt_unsc(p).p;
        Results_tracking(ww).ww(p).Qs_toTrack = Qs(p).p.allinterpfilt;
        Results_tracking(ww).ww(p).Ts_toTrack = ID(p).p.allinterp;
        Results_tracking(ww).ww(p).GRFs_toTrack = GRF(p).p.val.allinterp;
        Results_tracking(ww).ww(p).GRMs_toTrack = GRF(p).p.MorGF.allinterp;     
        Results_tracking(ww).ww(p).ParamsCM_opt = full(paramsCM_opt_nsc); 
        Results_tracking(ww).ww(p).bounds.params.upper = full(f_nscnp(...
            bounds(p).p.params.upper,scaling(p).p.params.v,...
            scaling(p).p.params.r));
        Results_tracking(ww).ww(p).bounds.params.lower = full(f_nscnp(...
            bounds(p).p.params.lower,scaling(p).p.params.v,...
            scaling(p).p.params.r));        
        Results_tracking(ww).ww(p).ParamsCM_gen = full(f_nscnp(...
            guess(1).p.params,scaling(1).p.params.v,scaling(1).p.params.r));
        Results_tracking(ww).ww(p).dev_cm = dev_cm;
        Results_tracking(ww).ww(p).stats = stats;
        Results_tracking(ww).ww(p).colheaders.joints = joints;
        Results_tracking(ww).ww(p).colheaders.GRF = {'fore_aft_r','vertical_r',...
            'lateral_r','fore_aft_l','vertical_l','lateral_l'};
        for i = 1:NMuscle/2
                Results_tracking(ww).ww(p).colheaders.muscles{i} = ...
                    [muscleNames{i}(1:end-2),'_l'];
                Results_tracking(ww).ww(p).colheaders.muscles{i+NMuscle/2} = ...
                    [muscleNames{i}(1:end-2),'_r'];
        end
        Results_tracking(ww).ww(p).colheaders.paramsCM = {'loc_s1_r_x','loc_s1_r_z',...
            'loc_s2_r_x','loc_s2_r_z','loc_s3_r_x','loc_s3_r_z',...
            'loc_s4_r_x','loc_s4_r_z','loc_s5_r_x','loc_s5_r_z',...
            'loc_s6_r_x','loc_s6_r_z','radius_s1','radius_s2','radius_s3',...
            'radius_s4','radius_s5','radius_s6'};
    end
    % Save data
    save([pathresults,'/',namescript,'/Results_tracking.mat'],...
        'Results_tracking');
    end
end
end
