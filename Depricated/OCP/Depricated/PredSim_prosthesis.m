%%  Three-dimensional muscle-driven predictive simulations of human gaits
%   with a passive transtibial prosthesis
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

% settings describes the parameters used in the optimal control problems.
% settings(1): average speed: max 2 digits.
% settings(2): NLP error tolerance: 1*10^(-settings(2)).
% settings(3): number of mesh intervals
% settings(4): weight metabolic energy rate
% settings(5): weight joint accelerations
% settings(6): weight arm excitations
% settings(7): weight passive torques
% settings(8): weight muscle activations
% settings(9): power metabolic energy rate
% settings(10): initial guess identifier: 1 quasi-random, 2 data-informed
% settings(11): contact model identifier: 1 generic, 2 subject-specific
% settings(12): initial guess mode identifier for data-informed guesses
%   1 data from average walking motion, 2 data from average running motion 
%   3 partial solution from motion specified in settings(13)
%   4 full solution from motion specified in settings(13)
%   if 1-3 then settings(10) should be 2
%   if 4 then settings(10) does not matter
% settings(13): initial guess case identifier: row number in "settings" of
% the motion used as initial guess when settings(12)=3 || settings(12)=4
% settings(14): weakness hip muscles, 0 no weakness, 50 weakened by 50%,
% 75 weakened by 75%, 90 weakened by 90%,
% settings(15): maximal contraction velocity identifier
%   0 generic values, 1 values*2
% settings(16): weakness ankle plantarflexors, 0 no weakness, 50 weakened 
% by 50%, 75 weakened by 75%, 90 weakened by 90%,
% settings(17): metabolic energy model identifier: 0 Bhargava et al. (2004)
% 1 Umberger et al. (2003), 2 Umberger (2010), 3 Uchida et al. (2016), 
% 4 Umberger (2010) treating muscle lengthening heat rate as 
% Umberger et al. (2003), 5 Umberger (2010) treating negative mechanical 
% work as Umberger et al. (2003)
% settings(18): torsional stiffness (N m rad-1)
settings = [        
    % Data-informed (walking) initial guess
    1.33, 4, 100, 500, 50000, 1000000, 1000, 2000, 2, 2, 1, 1, 0, 0, 0, 0, 0, 800;   % 1
    % Quasi-random initial guess
    1.33, 4, 100, 500, 50000, 1000000, 1000, 2000, 2, 1, 1, 0, 0, 0, 0, 0, 0, 800;   % 2
    ];

%% Select settings
for www = 1:length(idx_ww)
    
%% Set parameters based on settings
ww = idx_ww(www);
%% Variable parameters
v_tgt       = settings(ww,1);    % target velocity
tol_ipopt   = settings(ww,2);    % tolerance (means 1e-(tol_ipopt))
N           = settings(ww,3);    % number of mesh intervals
W.E         = settings(ww,4);    % weight metabolic energy
W.Ak        = settings(ww,5);    % weight joint accelerations
W.ArmE      = settings(ww,6);    % weight arm excitations
W.passMom   = settings(ww,7);    % weight passive moments
W.A         = settings(ww,8);    % weight muscle activations
exp_E       = settings(ww,9);    % power metabolic energy
IGsel       = settings(ww,10);   % initial guess selection
cm          = settings(ww,11);   % contact model
IGm         = settings(ww,12);   % initial guess mode
IGcase      = settings(ww,13);   % initial guess case
h_weak      = settings(ww,14);   % h_weakness hip actuators
vMax_s      = settings(ww,15);   % maximal contraction velocity
pf_weak     = settings(ww,16);   % weakness ankle plantaflexors
mE          = settings(ww,17);   % metabolic energy model identifier
kstiff      = settings(ww,18);   % prosthesis stiffness
% Fixed parameters
W.u = 0.001;
% The filename used to save the results depends on the settings 
v_tgt_id = round(v_tgt,2);
savename = ['_c',num2str(ww),'_v',num2str(v_tgt_id*100),...
    '_T',num2str(tol_ipopt),'_N',num2str(N),'_E',num2str(W.E),...
    '_Ak',num2str(W.Ak),'_AE',num2str(W.ArmE),'_P',num2str(W.passMom),...
    '_A',num2str(W.A),'_eE',num2str(exp_E),'_G',num2str(IGsel),...
    '_M',num2str(cm),'_Gm',num2str(IGm),...
    '_W',num2str(h_weak),'_vM',num2str(vMax_s),...
    '_pW',num2str(pf_weak),'_mE',num2str(mE),'_k',num2str(kstiff)];

%% Load external functions
% The external function performs inverse dynamics through the
% OpenSim/Simbody C++ API. This external function is compiled as a dll from
% which we create a Function instance using CasADi in MATLAB.
pathmain = pwd;
% We use different external functions, since we also want to access some 
% parameters of the model in a post-processing phase.
[pathRepo,~,~] = fileparts(pathmain);
pathExternalFunctions = [pathRepo,'/ExternalFunctions'];                 
% Loading external functions. 
cd(pathExternalFunctions);
setup.derivatives =  'AD'; % Algorithmic differentiation
if ispc    
    switch setup.derivatives
        case {'AD'}
            F = external('F','PredSimProsthesis.dll');   
            if analyseResults
                F1 = external('F','PredSimProsthesis_pp.dll');
            end
    end
elseif ismac
    switch setup.derivatives
        case {'AD'}
            F = external('F','PredSimProsthesis.dylib');   
            if analyseResults
                F1 = external('F','PredSimProsthesis_pp.dylib');
            end
    end   
else
    disp('Platform not supported')
end
cd(pathmain);
% This is an example of how to call an external function with some
% numerical values.
% vec1 = -ones(87,1);
% res1 = full(F(vec1));
% res2 = full(F1(vec1));

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
residuals_noarmsi   = jointi.pelvis.tilt:jointi.trunk.rot; % all but arms
% Number of degrees of freedom for later use
nq.all      = length(residualsi); % all 
nq.abs      = length(ground_pelvisi); % ground-pelvis
nq.trunk    = length(trunki); % trunk
nq.arms     = length(armsi); % arms
nq.leg      = 9; % #joints needed for polynomials
% Second, origins bodies. 
% Calcaneus
calcOr.r    = 30:31;
calcOr.l    = 32:33;
calcOr.all  = [calcOr.r,calcOr.l];
NcalcOr     = length(calcOr.all);
% Femurs
femurOr.r   = 34:35;
femurOr.l   = 36:37;
femurOr.all = [femurOr.r,femurOr.l];
NfemurOr    = length(femurOr.all);
% Hands
handOr.r    = 38:39;
handOr.l    = 40:41;
handOr.all  = [handOr.r,handOr.l];
NhandOr     = length(handOr.all);
% Tibias
tibiaOr.r   = 42:43;
tibiaOr.l   = 44:45;
tibiaOr.all = [tibiaOr.r,tibiaOr.l];
NtibiaOr    = length(tibiaOr.all);
% External function: F1 (post-processing purpose only)
% Ground reaction forces (GRFs)
GRFi.r      = 30:32;
GRFi.l      = 33:35;
GRFi.all    = [GRFi.r,GRFi.l];
NGRF        = length(GRFi.all);
% Origins calcaneus (3D)
calcOrall.r     = 36:38;
calcOrall.l     = 39:41;
calcOrall.all   = [calcOrall.r,calcOrall.l];
NcalcOrall      = length(calcOrall.all);

%% Model info
body_mass = 8.84259166189724+2*6.98382288222561+2.78372323906632+...
    1.809420105393108+0.0750835667988218+0.048804318419234+...
    0.938544584985273+0.610053980240428+0.162631005686248+...
    0.105710153696061+25.7060604306454+2*1.52611854532613+...
    2*0.456132668302843+2*0.456132668302843+2*0.34350731810461;
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
% Muscles replaced by prosthesis    
muscleNames_pro = {'med_gas_r','lat_gas_r','soleus_r','tib_post_r',...
    'flex_dig_r','flex_hal_r','tib_ant_r','per_brev_r','per_long_r',...
    'per_tert_r','ext_dig_r','ext_hal_r'};  
% Muscles from right side (without muscles replaced by prosthesis)
muscleNames_act = {'glut_med1_r','glut_med2_r','glut_med3_r',...
        'glut_min1_r','glut_min2_r','glut_min3_r','semimem_r',...
        'semiten_r','bifemlh_r','bifemsh_r','sar_r','add_long_r',...
        'add_brev_r','add_mag1_r','add_mag2_r','add_mag3_r','tfl_r',...
        'pect_r','grac_r','glut_max1_r','glut_max2_r','glut_max3_r',......
        'iliacus_r','psoas_r','quad_fem_r','gem_r','peri_r',...
        'rect_fem_r','vas_med_r','vas_int_r','vas_lat_r',...
        'ercspn_r','intobl_r','extobl_r'};
% Muscle indices for later use    
pathmusclemodel = [pathRepo,'/MuscleModel'];
addpath(genpath(pathmusclemodel));    
% (1:end-3), since we do not want to count twice the back muscles
musi = MuscleIndices(muscleNames(1:end-3));
musi_pro = MuscleIndices(muscleNames_act);
% Total number of muscles
NMuscle = length(muscleNames(1:end-3))*2;
NMuscle_pro = length(muscleNames_pro);
NMuscle_act = NMuscle-NMuscle_pro;
% Muscle-tendon parameters. Row 1: maximal isometric forces; Row 2: optimal
% fiber lengths; Row 3: tendon slack lengths; Row 4: optimal pennation 
% angles; Row 5: maximal contraction velocities
load([pathmusclemodel,'/MTparameters_',subject,'.mat']);
MTparameters_m = [MTparameters(:,musi),MTparameters(:,musi_pro)];
% Indices of the muscles actuating the different joints for later use
pathpolynomial = [pathRepo,'/Polynomials'];
addpath(genpath(pathpolynomial));
tl = load([pathpolynomial,'/muscle_spanning_joint_INFO_',subject,'.mat']);
[~,mai] = MomentArmIndices(muscleNames(1:end-3),...
    tl.muscle_spanning_joint_INFO(1:end-3,:));
% Prosthesis. For the right knee, no moment arms from the gastrocs
% Hard coded: indices 12 and 13 in mai(4).mus.l correspond to med gas (32)
% and lat gas (33). We need to take the indices from the left side since
% it is used with the polynomials that are calculated leg by leg
mai(4).mus.r_pro = mai(4).mus.l;
mai(4).mus.r_pro([12,13]) = [];
mai(4).mus.r([12,13]) = [];
% For the back, we also need to adjust the order
mai(7).mus.r = [78,79,80];
mai(8).mus.r = [78,79,80];
mai(9).mus.r = [78,79,80];

% By default, the tendon stiffness is 35 and the shift is 0.
aTendon = 35*ones(NMuscle_act,1);
shift = zeros(NMuscle_act,1);

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

%% Metabolic energy model parameters
% We extract the specific tensions and slow twitch rations.
pathMetabolicEnergy = [pathRepo,'/MetabolicEnergy'];
addpath(genpath(pathMetabolicEnergy));
% (1:end-3), since we do not want to count twice the back muscles
tension = getSpecificTensions(muscleNames(1:end-3)); 
tensions = [tension;tension(musi_pro)];
% (1:end-3), since we do not want to count twice the back muscles
pctst = getSlowTwitchRatios(muscleNames(1:end-3)); 
pctsts = [pctst;pctst(musi_pro)];

%% CasADi functions
% We create several CasADi functions for later use
pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
addpath(genpath(pathCasADiFunctions));
% We load some variables for the polynomial approximations
load([pathpolynomial,'/muscle_spanning_joint_INFO_',subject,'.mat']);
load([pathpolynomial,'/MuscleInfo_',subject,'.mat']);
% For the polynomials, we want all independent muscles. So we do not need 
% the muscles from both legs, since we assume bilateral symmetry, but want
% all muscles from the back (indices 47:49).
musi_pol = [musi,47,48,49];
NMuscle_pol = NMuscle/2+3;
CasADiFunctions_prosthesis

%% Passive joint torques
% We extract the parameters for the passive torques of the lower limbs and
% the trunk
pathPassiveMoments = [pathRepo,'/PassiveMoments'];
addpath(genpath(pathPassiveMoments));
PassiveMomentsData

%% Experimental data
% We extract experimental data to set bounds and initial guesses if needed
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
% Extract joint positions from average walking motion
motion_walk         = 'walking';
nametrial_walk.id   = ['average_',motion_walk,'_FGC']; 
nametrial_walk.IK	= ['IK_',nametrial_walk.id];
pathIK_walk         = [pathData,'/IK/',nametrial_walk.IK,'.mat'];
Qs_walk             = getIK(pathIK_walk,joints);  

%% Bounds
pathBounds = [pathRepo,'/Bounds'];
addpath(genpath(pathBounds));
[bounds,scaling] = getBounds_prosthesis(Qs_walk,NMuscle_act,nq,jointi);

%% Initial guess
% The initial guess depends on the settings
pathIG = [pathRepo,'/IG'];
addpath(genpath(pathIG));
if IGsel == 1  % Quasi-random initial guess  
    guess = getGuess_QR_prosthesis(N,nq,NMuscle_act,scaling,v_tgt,jointi);
elseif IGsel == 2 % Data-informed initial guess
    if IGm == 1 % Data from average walking motion
        time_IC = [Qs_walk.time(1),Qs_walk.time(end)];
        guess = getGuess_DI_prosthesis(Qs_walk,nq,N,time_IC,...
            NMuscle_act,jointi,scaling,v_tgt);   
    end 
end
% This allows visualizing the initial guess and the bounds
if checkBoundsIG
    pathPlots = [pathRepo,'/Plots'];
    addpath(genpath(pathPlots));
    plot_BoundsVSInitialGuess_prosthesis
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
    % Final time
    tf              = MX.sym('tf',1);
    w               = [w {tf}];
    lbw             = [lbw; bounds.tf.lower];
    ubw             = [ubw; bounds.tf.upper];
    w0              = [w0;  guess.tf];
    % Define states at first mesh point
    % Muscle activations
    a0              = MX.sym('a0',NMuscle_act);
    w               = [w {a0}];
    lbw             = [lbw; bounds.a.lower'];
    ubw             = [ubw; bounds.a.upper'];
    w0              = [w0;  guess.a(1,:)'];
    % Muscle-tendon forces
    FTtilde0        = MX.sym('FTtilde0',NMuscle_act);
    w               = [w {FTtilde0}];
    lbw             = [lbw; bounds.FTtilde.lower'];
    ubw             = [ubw; bounds.FTtilde.upper'];
    w0              = [w0;  guess.FTtilde(1,:)'];
    % Qs and Qdots 
    X0              = MX.sym('X0',2*nq.all);
    w               = [w {X0}];    
    lbw             = [lbw; bounds.QsQdots_0.lower'];
    ubw             = [ubw; bounds.QsQdots_0.upper'];    
    w0              = [w0;  guess.QsQdots(1,:)'];
    % Arm activations
    a_a0            = MX.sym('a_a0',nq.arms);
    w               = [w {a_a0}];
    lbw             = [lbw; bounds.a_a.lower'];
    ubw             = [ubw; bounds.a_a.upper'];
    w0              = [w0;  guess.a_a(1,:)'];
    % We pre-allocate some of the states so that we can provide an
    % expression for the distance traveled  
    for k=0:N
        Xk{k+1,1} = MX.sym(['X_' num2str(k+1)], 2*nq.all);
    end 
    % "Lift" initial conditions
    ak          = a0;
    FTtildek    = FTtilde0;
    Xk{1,1}     = X0;
    a_ak        = a_a0; 
    % Provide expression for the distance traveled
    pelvis_tx0 = Xk{1,1}(2*jointi.pelvis.tx-1,1).*...
        scaling.QsQdots(2*jointi.pelvis.tx-1); % initial position pelvis_tx     
    pelvis_txf = Xk{N+1,1}(2*jointi.pelvis.tx-1,1).*...
        scaling.QsQdots(2*jointi.pelvis.tx-1); % final position pelvis_tx    
    dist_trav_tot = pelvis_txf-pelvis_tx0;% distance traveled  
    % Time step
    h = tf/N; 
    % Loop over mesh points
    for k=0:N-1
        % Define controls at mesh point (piecewise-constant in interval) 
        % Time derivative of muscle activations (states)
        vAk         = MX.sym(['vA_' num2str(k)], NMuscle_act);
        w           = [w {vAk}];
        lbw         = [lbw; bounds.vA.lower'];
        ubw         = [ubw; bounds.vA.upper'];
        w0          = [w0; guess.vA(k+1,:)'];
        % Time derivative of muscle-tendon forces (states)
        dFTtildek   = MX.sym(['dFTtilde_' num2str(k)], NMuscle_act);
        w           = [w {dFTtildek}];
        lbw         = [lbw; bounds.dFTtilde.lower'];
        ubw         = [ubw; bounds.dFTtilde.upper'];
        w0          = [w0; guess.dFTtilde(k+1,:)'];  
        % Time derivative of Qdots (states) 
        Ak          = MX.sym(['A_' num2str(k)], nq.all);
        w           = [w {Ak}];
        lbw         = [lbw; bounds.Qdotdots.lower'];
        ubw         = [ubw; bounds.Qdotdots.upper'];
        w0          = [w0; guess.Qdotdots(k+1,:)'];   
        % Arm excitations
        e_ak        = MX.sym(['e_a_' num2str(k)], nq.arms);
        w           = [w {e_ak}];
        lbw         = [lbw; bounds.e_a.lower'];
        ubw         = [ubw; bounds.e_a.upper'];
        w0          = [w0; guess.e_a(k+1,:)'];
        % Define states at collocation points    
        % Muscle activations
        akj = {};
        for j=1:d
            akj{j}=MX.sym(['a_' num2str(k) '_' num2str(j)],NMuscle_act);
            w     = {w{:}, akj{j}};
            lbw   = [lbw; bounds.a.lower'];
            ubw   = [ubw; bounds.a.upper'];
            w0    = [w0;  guess.a(k+1,:)'];
        end   
        % Muscle-tendon forces
        FTtildekj = {};
        for j=1:d
            FTtildekj{j} = ...
                MX.sym(['FTtilde_' num2str(k) '_' num2str(j)],NMuscle_act);
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
            a_akj{j}= MX.sym(['a_a_' num2str(k) '_' num2str(j)], nq.arms);
            w       = {w{:}, a_akj{j}};
            lbw     = [lbw; bounds.a_a.lower'];
            ubw     = [ubw; bounds.a_a.upper'];
            w0      = [w0;  guess.a_a(k+1,:)'];
        end   
        % Unscale variables for later use
        Xk_nsc          = Xk{k+1,1}.*scaling.QsQdots';
        FTtildek_nsc    = FTtildek.*(scaling.FTtilde');
        Ak_nsc          = Ak.*scaling.Qdotdots';
        for j=1:d
            Xkj_nsc{j} = Xkj{j}.*scaling.QsQdots';
            FTtildekj_nsc{j} = FTtildekj{j}.*scaling.FTtilde';
        end   
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
        MA.knee.r       =  MA_r(mai(4).mus.r_pro',4);
        % Both legs
        % In MuscleInfo, we first have the right back muscles (44:46) and 
        % then the left back muscles (47:49). Here we re-organize so that
        % we have first the left muscles and then the right muscles.
        % Prosthesis: we remove the prosthesis muscles for the right leg
        lMTk_lr     = [lMTk_l([1:43,47:49],1);lMTk_r(musi_pro,1)];
        vMTk_lr     = [vMTk_l([1:43,47:49],1);vMTk_r(musi_pro,1)];   
        % Get muscle-tendon forces and derive Hill-equilibrium
        [Hilldiffk,FTk,Fcek,Fpassk,Fisok,vMmaxk,massMk] = ...
            f_forceEquilibrium_FtildeState_all_tendon(...
                ak,FTtildek.*scaling.FTtilde',...
                dFTtildek.*scaling.dFTtilde,lMTk_lr,vMTk_lr,tensions,...
                aTendon,shift);   
        % Get metabolic energy rate if in the cost function  
        if W.E ~= 0   
            % Get muscle fiber lengths
            [~,lMtildek] = f_FiberLength_TendonForce_tendon(...
                FTtildek.*scaling.FTtilde',lMTk_lr,aTendon,shift); 
            % Get muscle fiber velocities
            [vMk,~] = f_FiberVelocity_TendonForce_tendon(...
                FTtildek.*scaling.FTtilde',...
                dFTtildek.*scaling.dFTtilde,lMTk_lr,vMTk_lr,aTendon,shift);
            % Bhargava et al. (2004) 
            [e_tot,~,~,~,~,~] = fgetMetabolicEnergySmooth2004all(...
                ak,ak,lMtildek,vMk,Fcek,Fpassk,massMk,...
                pctsts,Fisok,MTparameters_m(1,:)',body_mass,10);          
        end
        % Get passive joint torques
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
        Tau_passk.subt.l       = f_PassiveMoments(k_pass.subt,...
            theta.pass.subt,Xk_nsc(jointi.subt.l*2-1,1),...
            Xk_nsc(jointi.subt.l*2,1));    
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
            Tau_passk.subt.l,...
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
            lbg     = [lbg; zeros(NMuscle_act,1)];
            ubg     = [ubg; zeros(NMuscle_act,1)]; 
            % Contraction dynamics (implicit formulation)              
            g       = {g{:}, (h*dFTtildek.*scaling.dFTtilde - ...
                FTtildep_nsc)./(scaling.FTtilde')};
            lbg     = [lbg; zeros(NMuscle_act,1)];
            ubg     = [ubg; zeros(NMuscle_act,1)];
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
            if W.E ~= 0
            J = J + 1/(dist_trav_tot)*(...
                W.E*B(j+1)      *(f_JNMactexp(e_tot,exp_E))/body_mass*h+...                
                W.A*B(j+1)      *(f_JNMact(akj{j}))*h + ...
                W.ArmE*B(j+1)   *(f_J8(e_ak))*h +... 
                W.Ak*B(j+1)     *(f_J21(Ak(residuals_noarmsi,1)))*h + ... 
                W.passMom*B(j+1)*(f_J13(Tau_passk_all))*h + ...
                W.u*B(j+1)      *(f_JNMact(vAk))*h + ...
                W.u*B(j+1)      *(f_JNMact(dFTtildek))*h + ...                
                W.u*B(j+1)      *(f_J8(Ak(armsi,1)))*h);
            else
            J = J + 1/(dist_trav_tot)*(...
                W.A*B(j+1)      *(f_JNMact(akj{j}))*h + ...
                W.ArmE*B(j+1)   *(f_J8(e_ak))*h +... 
                W.Ak*B(j+1)     *(f_J21(Ak(residuals_noarmsi,1)))*h + ... 
                W.passMom*B(j+1)*(f_J13(Tau_passk_all))*h + ...
                W.u*B(j+1)      *(f_JNMact(vAk))*h + ...
                W.u*B(j+1)      *(f_JNMact(dFTtildek))*h + ...                
                W.u*B(j+1)      *(f_J8(Ak(armsi,1)))*h);
            end
        end            
        % Call external function
        [Tk] = F([Xk_nsc;Ak_nsc]);            
        % Add path constraints
        % Null pelvis residuals
        g               = {g{:},Tk(ground_pelvisi,1)};
        lbg             = [lbg; zeros(nq.abs,1)];
        ubg             = [ubg; zeros(nq.abs,1)];    
        % Muscle-driven joint torques for the lower limbs and the trunk
        % Hip flexion, left
        Ft_hip_flex_l   = FTk(mai(1).mus.l',1);
        T_hip_flex_l    = f_T27(MA.hip_flex.l,Ft_hip_flex_l);
        g               = {g{:},Tk(jointi.hip_flex.l,1)-(T_hip_flex_l + ...
            Tau_passk.hip.flex.l - 0.1*Xk_nsc(jointi.hip_flex.l*2,1))};        
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];    
        % Hip flexion, right
        Ft_hip_flex_r   = FTk(mai(1).mus.r',1);
        T_hip_flex_r    = f_T27(MA.hip_flex.r,Ft_hip_flex_r);
        g               = {g{:},Tk(jointi.hip_flex.r,1)-(T_hip_flex_r + ...
            Tau_passk.hip.flex.r - 0.1*Xk_nsc(jointi.hip_flex.r*2,1))};        
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];   
        % Hip adduction, left
        Ft_hip_add_l    = FTk(mai(2).mus.l',1);
        T_hip_add_l     = f_T27(MA.hip_add.l,Ft_hip_add_l);
        g               = {g{:},Tk(jointi.hip_add.l,1)-(T_hip_add_l + ...
            Tau_passk.hip.add.l - 0.1*Xk_nsc(jointi.hip_add.l*2,1))};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];    
        % Hip adduction, right
        Ft_hip_add_r    = FTk(mai(2).mus.r',1);
        T_hip_add_r     = f_T27(MA.hip_add.r,Ft_hip_add_r);
        g               = {g{:},Tk(jointi.hip_add.r,1)-(T_hip_add_r + ...
            Tau_passk.hip.add.r - 0.1*Xk_nsc(jointi.hip_add.r*2,1))};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];  
        % Hip rotation, left
        Ft_hip_rot_l    = FTk(mai(3).mus.l',1);
        T_hip_rot_l     = f_T27(MA.hip_rot.l,Ft_hip_rot_l);
        g               = {g{:},Tk(jointi.hip_rot.l,1)-(T_hip_rot_l + ...
            Tau_passk.hip.rot.l - 0.1*Xk_nsc(jointi.hip_rot.l*2,1))};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];    
        % Hip rotation, right
        Ft_hip_rot_r    = FTk(mai(3).mus.r',1);
        T_hip_rot_r     = f_T27(MA.hip_rot.r,Ft_hip_rot_r);
        g               = {g{:},Tk(jointi.hip_rot.r,1)-(T_hip_rot_r + ...
            Tau_passk.hip.rot.r - 0.1*Xk_nsc(jointi.hip_rot.r*2,1))};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];   
        % Knee, left
        Ft_knee_l       = FTk(mai(4).mus.l',1);
        T_knee_l        = f_T13(MA.knee.l,Ft_knee_l);
        g               = {g{:},Tk(jointi.knee.l,1)-(T_knee_l + ...
            Tau_passk.knee.l - 0.1*Xk_nsc(jointi.knee.l*2,1))};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];    
        % Knee, right
        Ft_knee_r       = FTk(mai(4).mus.r',1);
        T_knee_r        = f_T11(MA.knee.r,Ft_knee_r);
        g               = {g{:},Tk(jointi.knee.r,1)-(T_knee_r + ...
            Tau_passk.knee.r - 0.1*Xk_nsc(jointi.knee.r*2,1))};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];    
        % Ankle, left
        Ft_ankle_l      = FTk(mai(5).mus.l',1);
        T_ankle_l       = f_T12(MA.ankle.l,Ft_ankle_l);
        g               = {g{:},Tk(jointi.ankle.l,1)-(T_ankle_l + ...
            Tau_passk.ankle.l - 0.1*Xk_nsc(jointi.ankle.l*2,1))};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];    
        % Prosthesis: passive ankle stiffness
        g  = {g{:},Tk(jointi.ankle.r,1) - ...
            (-kstiff*Xk_nsc(jointi.ankle.r*2-1,1))};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];        
        % Subtalar, left
        Ft_subt_l       = FTk(mai(6).mus.l',1);
        T_subt_l        = f_T12(MA.subt.l,Ft_subt_l);
        g               = {g{:},(Tk(jointi.subt.l,1)-(T_subt_l + ...
            Tau_passk.subt.l - 0.1*Xk_nsc(jointi.subt.l*2,1)))};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];   
        % Prosthesis: passive subtalar stiffness
        g = {g{:},(Tk(jointi.subt.r,1) - ...
            (-kstiff*Xk_nsc(jointi.subt.r*2-1,1)))};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];
        % Lumbar extension
        Ft_trunk_ext    = FTk([mai(7).mus.l,mai(7).mus.r]',1);
        T_trunk_ext     = f_T6(MA.trunk_ext,Ft_trunk_ext);
        g               = {g{:},Tk(jointi.trunk.ext,1)-(T_trunk_ext + ...
            Tau_passk.trunk.ext - 0.1*Xk_nsc(jointi.trunk.ext*2,1))};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];
        % Lumbar bending
        Ft_trunk_ben    = FTk([mai(8).mus.l,mai(8).mus.r]',1);
        T_trunk_ben     = f_T6(MA.trunk_ben,Ft_trunk_ben);
        g               = {g{:},Tk(jointi.trunk.ben,1)-(T_trunk_ben + ...
            Tau_passk.trunk.ben - 0.1*Xk_nsc(jointi.trunk.ben*2,1))};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];
        % Lumbar rotation
        Ft_trunk_rot    = FTk([mai(9).mus.l,mai(9).mus.r]',1);
        T_trunk_rot     = f_T6(MA.trunk_rot,Ft_trunk_rot);
        g               = {g{:},Tk(jointi.trunk.rot,1)-(T_trunk_rot + ...
            Tau_passk.trunk.rot - 0.1*Xk_nsc(jointi.trunk.rot*2,1))};
        lbg             = [lbg; 0];
        ubg             = [ubg; 0];
        % Torque-driven joint torques for the arms
        % Arms
        g               = {g{:},Tk(armsi,1)/scaling.ArmTau - a_ak + ...
            0.1/scaling.ArmTau*Xk_nsc(armsi*2,1)};
        lbg             = [lbg; zeros(nq.arms,1)];
        ubg             = [ubg; zeros(nq.arms,1)];
        % Activation dynamics (implicit formulation)
        tact = 0.015;
        tdeact = 0.06;
        act1 = vAk*scaling.vA + ak./(ones(size(ak,1),1)*tdeact);
        act2 = vAk*scaling.vA + ak./(ones(size(ak,1),1)*tact);
        % act1
        g               = {g{:},act1};
        lbg             = [lbg; zeros(NMuscle_act,1)];
        ubg             = [ubg; inf*ones(NMuscle_act,1)]; 
        % act2
        g               = {g{:},act2};
        lbg             = [lbg; -inf*ones(NMuscle_act,1)];
        ubg             = [ubg; ones(NMuscle_act,1)./...
            (ones(NMuscle_act,1)*tact)];        
        % Contraction dynamics (implicit formulation)
        g               = {g{:},Hilldiffk};
        lbg             = [lbg; zeros(NMuscle_act,1)];
        ubg             = [ubg; zeros(NMuscle_act,1)];  
        % Constraints to prevent parts of the skeleton to penetrate each
        % other.
        % Origins calcaneus (transv plane) at minimum 9 cm from each other.    
        g               = {g{:},f_Jnn2(Tk(calcOr.r,1) - Tk(calcOr.l,1))};
        lbg             = [lbg; 0.0081];
        ubg             = [ubg; 4];   
        % Constraint to prevent the arms to penetrate the skeleton
        % Squared Euclidian distance (x-z plane) between origins femur 
        % and hand larger than 0.18^2 = 0.0324  
        g               = {g{:},f_Jnn2(Tk(femurOr.r,1) - Tk(handOr.r,1))};
        lbg             = [lbg; 0.0324];
        ubg             = [ubg; 4]; 
        % Constraint to prevent the arms to penetrate the skeleton       
        % Origins femurs and ipsilateral hands (transv plane) at minimum 
        % 18 cm from each other.    
        g               = {g{:},f_Jnn2(Tk(femurOr.l,1) - Tk(handOr.l,1))};
        lbg             = [lbg; 0.0324];
        ubg             = [ubg; 4]; 
        % Origins tibia (transv plane) at minimum 11 cm from each other.    
        g               = {g{:},f_Jnn2(Tk(tibiaOr.r,1) - Tk(tibiaOr.l,1))};
        lbg             = [lbg; 0.0121];
        ubg             = [ubg; 4];   
        % New NLP variables for states at end of interval
        if k ~= N-1
            % Muscle activations
            ak          = MX.sym(['a_' num2str(k+1)], NMuscle_act);
            w           = {w{:}, ak};
            lbw         = [lbw; bounds.a.lower'];
            ubw         = [ubw; bounds.a.upper'];
            w0          = [w0;  guess.a(k+2,:)'];
            % Muscle-tendon forces
            FTtildek    = MX.sym(['FTtilde_' num2str(k+1)], NMuscle_act);
            w           = {w{:}, FTtildek};
            lbw         = [lbw; bounds.FTtilde.lower'];
            ubw         = [ubw; bounds.FTtilde.upper'];
            w0          = [w0;  guess.FTtilde(k+2,:)'];    
            % Qs and Qdots
            w           = {w{:}, Xk{k+2,1}};
            lbw         = [lbw; bounds.QsQdots.lower'];
            ubw         = [ubw; bounds.QsQdots.upper']; 
            w0          = [w0;  guess.QsQdots(k+2,:)'];
            % Arm activations
            a_ak        = MX.sym(['a_a_' num2str(k+1)], nq.arms);
            w           = {w{:}, a_ak};
            lbw         = [lbw; bounds.a_a.lower'];
            ubw         = [ubw; bounds.a_a.upper'];
            w0          = [w0;  guess.a_a(k+2,:)'];
        else
            % Muscle activations
            ak          = MX.sym(['a_' num2str(k+1)], NMuscle_act);
            w           = {w{:}, ak};
            lbw         = [lbw; bounds.a.lower'];
            ubw         = [ubw; bounds.a.upper'];            
            w0          = [w0;  guess.a(1,:)'];
            % Muscle-tendon forces
            FTtildek    = MX.sym(['FTtilde_' num2str(k+1)], NMuscle_act);
            w           = {w{:}, FTtildek};
            lbw         = [lbw; bounds.FTtilde.lower'];
            ubw         = [ubw; bounds.FTtilde.upper'];
            w0          = [w0;  guess.FTtilde(1,:)'];    
            % Qs and Qdots
            w            = {w{:}, Xk{k+2,1}};
            lbw          = [lbw; bounds.QsQdots.lower'];
            ubw          = [ubw; bounds.QsQdots.upper'];
            end_X        = guess.QsQdots(1,:);              
            dx = guess.QsQdots(end,2*jointi.pelvis.tx-1) - ...
                guess.QsQdots(end-1,2*jointi.pelvis.tx-1);            
            end_X(2*jointi.pelvis.tx-1) = ...
                guess.QsQdots(end,2*jointi.pelvis.tx-1) + dx;            
            w0           = [w0;  end_X'];   
            % Arm activations
            a_ak         = MX.sym(['a_a_' num2str(k+1)], nq.arms);
            w            = {w{:}, a_ak};
            lbw          = [lbw; bounds.a_a.lower'];
            ubw          = [ubw; bounds.a_a.upper'];            
            w0           = [w0;  guess.a_a(1,:)'];
        end
        % Rescale variables to impose equality constraints
        Xk_end = (Xk_nsc_end)./scaling.QsQdots';
        FTtildek_end = (FTtildek_nsc_end)./scaling.FTtilde';
        % Add equality constraints (next interval starts with end values of 
        % states from previous interval)
        g   = {g{:}, Xk_end-Xk{k+2,1}, FTtildek_end-FTtildek, ...
            ak_end-ak, a_ak_end-a_ak};
        lbg = [lbg;zeros(2*nq.all+ NMuscle_act + NMuscle_act + nq.arms,1)];
        ubg = [ubg;zeros(2*nq.all+ NMuscle_act + NMuscle_act + nq.arms,1)]; 
    end  
    % Additional path constraints
    % Periodicity of the states
    % Qs and Qdots
    Qsper = [jointi.pelvis.tilt:2*jointi.pelvis.rot,...
        2*jointi.pelvis.tx:2*jointi.elb.r]';    
    g   = {g{:}, Xk_end(Qsper)-X0(Qsper,1)};
    lbg = [lbg; zeros(length(Qsper),1)];
    ubg = [ubg; zeros(length(Qsper),1)];       
    % Muscle activations
    g   = {g{:}, ak_end-a0};
    lbg = [lbg; zeros(NMuscle_act,1)];
    ubg = [ubg; zeros(NMuscle_act,1)];    
    % Muscle-tendon forces
    g   = {g{:}, FTtildek_end-FTtilde0};
    lbg = [lbg; zeros(NMuscle_act,1)];
    ubg = [ubg; zeros(NMuscle_act,1)];    
    % Arm activations
    g   = {g{:}, a_ak_end-a_a0};
    lbg = [lbg; zeros(nq.arms,1)];
    ubg = [ubg; zeros(nq.arms,1)];
    % Average speed
    vel_aver_tot = dist_trav_tot/tf; 
    g   = {g{:}, vel_aver_tot - v_tgt};
    lbg = [lbg; 0];
    ubg = [ubg; 0];
    
    % Assert bounds / IG
    % Lower bounds smaller than upper bounds    
    assert_bw = isempty(find(lbw <= ubw == 0,1)); % Design variables
    assert_bg = isempty(find(lbg <= ubg == 0,1)); % Constraints
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
    % Data-informed (full solution at closest speed) initial guess
    if IGm == 4
        load([pathresults,'/',namescript,'/w',savename_ig]);
        w0 = w_opt; clear w_opt; 
    end
    % Solve problem
    sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
        'lbg', lbg, 'ubg', ubg);    
    diary off
    % Extract results
    w_opt = full(sol.x);
    g_opt = full(sol.g);  
    % Save results and setup
    setup.tolerance.ipopt = tol_ipopt;
    setup.bounds = bounds;
    setup.scaling = scaling;
    setup.guess = guess;
    setup.lbw = lbw;
    setup.ubw = ubw;
    % Save results
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
    NControls = NMuscle_act+NMuscle_act+nq.all+nq.arms;
    NStates = NMuscle_act+NMuscle_act+2*nq.all+nq.arms;
    NParameters = 1;
    % In the loop
    Nwl = NControls+d*(NStates)+NStates;
    % In total
    Nw = NParameters+NStates+N*Nwl;
    % Before the variable corresponding to the first collocation point
    Nwm = NParameters+NStates+NControls;
    % Here we extract the results and re-organize them for analysis  
    % Static parameters
    tf_opt  = w_opt(1:NParameters);
    % Mesh points
    % Muscle activations and muscle-tendon forces
    a_opt = zeros(N+1,NMuscle_act);
    FTtilde_opt = zeros(N+1,NMuscle_act);
    for i = 1:NMuscle_act
        a_opt(:,i) = w_opt(NParameters+i:Nwl:Nw);
        FTtilde_opt(:,i) = w_opt(NParameters+NMuscle_act+i:Nwl:Nw);
    end
    % Qs and Qdots
    q_opt = zeros(N+1,nq.all);
    qdot_opt = zeros(N+1,nq.all);
    count = 0;
    for i = 1:2:2*nq.all
        count = count +1;
        q_opt(:,count) = ...
            w_opt(NParameters+NMuscle_act+NMuscle_act+i:Nwl:Nw);
        qdot_opt(:,count) = ...
            w_opt(NParameters+NMuscle_act+NMuscle_act+i+1:Nwl:Nw);
    end
    % Arm activations
    a_a_opt = zeros(N+1,nq.arms);
    for i = 1:nq.arms
        a_a_opt(:,i) = ...
            w_opt(NParameters+NMuscle_act+NMuscle_act+2*nq.all+i:Nwl:Nw);
    end    
    % Time derivative of muscle activations and muscle-tendon forces
    vA_opt = zeros(N,NMuscle_act);
    dFTtilde_opt = zeros(N,NMuscle_act);
    for i = 1:NMuscle_act
        vA_opt(:,i) = w_opt(NParameters+NStates+i:Nwl:Nw);
        dFTtilde_opt(:,i) =w_opt(NParameters+NStates+NMuscle_act+i:Nwl:Nw);
    end
    % Time derivative of joint velocities
    qdotdot_opt = zeros(N,nq.all);
    for i = 1:nq.all
        qdotdot_opt(:,i) = ...
            w_opt(NParameters+NStates+NMuscle_act+NMuscle_act+i:Nwl:Nw);
    end
    % Arm excitations
    e_a_opt = zeros(N,nq.arms);
    for i = 1:nq.arms
        e_a_opt(:,i) = w_opt(NParameters+NStates+NMuscle_act+...
            NMuscle_act+nq.all+i:Nwl:Nw); 
    end    
    % Collocation points
    % Muscle activations
    a_opt_ext=zeros(N*(d+1)+1,NMuscle_act);
    a_opt_ext(1:(d+1):end,:)= a_opt;
    for nmusi=1:NMuscle_act
        a_opt_ext(2:(d+1):end,nmusi) = w_opt(Nwm+nmusi:Nwl:Nw);
        a_opt_ext(3:(d+1):end,nmusi) = ...
            w_opt(Nwm+NMuscle_act+nmusi:Nwl:Nw);
        a_opt_ext(4:(d+1):end,nmusi) = ...
            w_opt(Nwm+NMuscle_act+NMuscle_act+nmusi:Nwl:Nw);
    end
    % Muscle activations at collocation points only
    a_opt_ext_col = zeros(N*d,NMuscle_act); 
    for nmusi=1:NMuscle_act
        a_opt_ext_col(1:d:end,nmusi) = w_opt(Nwm+nmusi:Nwl:Nw);
        a_opt_ext_col(2:d:end,nmusi) = w_opt(Nwm+NMuscle_act+nmusi:Nwl:Nw);
        a_opt_ext_col(3:d:end,nmusi) = ...
            w_opt(Nwm+NMuscle_act+NMuscle_act+nmusi:Nwl:Nw);   
    end    
    % Muscle-tendon forces
    FTtilde_opt_ext=zeros(N*(d+1)+1,NMuscle_act);
    FTtilde_opt_ext(1:(d+1):end,:)= FTtilde_opt;
    for nmusi=1:NMuscle_act
        FTtilde_opt_ext(2:(d+1):end,nmusi) = ...
            w_opt(Nwm+d*NMuscle_act+nmusi:Nwl:Nw);
        FTtilde_opt_ext(3:(d+1):end,nmusi) = ...
            w_opt(Nwm+d*NMuscle_act+NMuscle_act+nmusi:Nwl:Nw);
        FTtilde_opt_ext(4:(d+1):end,nmusi) = ...
            w_opt(Nwm+d*NMuscle_act+NMuscle_act+NMuscle_act+nmusi:Nwl:Nw);
    end
    % Qs and Qdots
    q_opt_ext=zeros(N*(d+1)+1,nq.all);
    q_opt_ext(1:(d+1):end,:)= q_opt;
    q_dot_opt_ext=zeros(N*(d+1)+1,nq.all);
    q_dot_opt_ext(1:(d+1):end,:)= qdot_opt;
    nqi_col = 1:2:2*nq.all;
    for nqi=1:nq.all
        nqi_q = nqi_col(nqi);
        q_opt_ext(2:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle_act+...
            d*NMuscle_act+nqi_q:Nwl:Nw);   
        q_opt_ext(3:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle_act+...
            d*NMuscle_act+2*nq.all+nqi_q:Nwl:Nw);  
        q_opt_ext(4:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle_act+...
            d*NMuscle_act+2*nq.all+2*nq.all+nqi_q:Nwl:Nw);  
        q_dot_opt_ext(2:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle_act+...
            d*NMuscle_act+nqi_q+1:Nwl:Nw);   
        q_dot_opt_ext(3:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle_act+...
            d*NMuscle_act+2*nq.all+nqi_q+1:Nwl:Nw);  
        q_dot_opt_ext(4:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle_act+...
            d*NMuscle_act+2*nq.all+2*nq.all+nqi_q+1:Nwl:Nw);
    end
    % Arm activations
    a_a_opt_ext=zeros(N*(d+1)+1,nq.arms);
    a_a_opt_ext(1:(d+1):end,:)= a_a_opt;
    for nmusi=1:nq.arms
        a_a_opt_ext(2:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle_act+...
            d*NMuscle_act+d*2*nq.all+nmusi:Nwl:Nw);
        a_a_opt_ext(3:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle_act+...
            d*NMuscle_act+d*2*nq.all+nq.arms+nmusi:Nwl:Nw);
        a_a_opt_ext(4:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle_act+...
            d*NMuscle_act+d*2*nq.all+nq.arms+nq.arms+nmusi:Nwl:Nw);
    end   
    
    %% Unscale results
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
        ones(1,NMuscle_act)*tdeact,ones(1,NMuscle_act)*tact);
    % Time derivative of muscle-tendon forces
    dFTtilde_opt_unsc = dFTtilde_opt.*repmat(...
        scaling.dFTtilde,size(dFTtilde_opt,1),size(dFTtilde_opt,2));
    % Arm excitations
    e_a_opt_unsc = e_a_opt.*repmat(scaling.e_a,size(e_a_opt,1),...
        size(e_a_opt,2));
    
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
 
    %% Joint torques and ground reaction forces at optimal solution
    Xk_Qs_Qdots_opt             = zeros(N,2*nq.all);
    Xk_Qs_Qdots_opt(:,1:2:end)  = q_opt_unsc.rad;
    Xk_Qs_Qdots_opt(:,2:2:end)  = qdot_opt_unsc.rad;
    Xk_Qdotdots_opt             = qdotdot_opt_unsc.rad;
    out_res_opt = zeros(N,nq.all+NGRF+NcalcOrall);
    for i = 1:N
        [res] = F1([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)']);
        out_res_opt(i,:) = full(res);    
    end
    GRF_opt_unsc = out_res_opt(:,GRFi.all);
    % assertArmTmax should be 0
    assertArmTmax = max(max(abs(out_res_opt(:,armsi)/scaling.ArmTau  -  ...
        (a_a_opt_unsc  + 0.1/scaling.ArmTau*Xk_Qs_Qdots_opt(:,armsi*2)))));
    
    %% Passive joint torques at optimal solution
    Tau_pass_opt_all = zeros(N,13);
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
        Tau_pass_opt.subt.l       = f_PassiveMoments(k_pass.subt,...
           theta.pass.subt,Xk_Qs_Qdots_opt(i,jointi.subt.l*2-1),...
           Xk_Qs_Qdots_opt(i,jointi.subt.l*2));
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
           Tau_pass_opt.subt.l,...
           Tau_pass_opt.trunk.ext,...
           Tau_pass_opt.trunk.ben,Tau_pass_opt.trunk.rot]);
    end    
    
    %% Assert average speed    
    dist_trav_opt = q_opt_ext(end,jointi.pelvis.tx)*...
        scaling.Qs(jointi.pelvis.tx) - q_opt_ext(1,jointi.pelvis.tx)*...
        scaling.Qs(jointi.pelvis.tx); % distance traveled
    time_elaps_opt = tf_opt; % time elapsed
    vel_aver_opt = dist_trav_opt/time_elaps_opt; 
    % assert_v_tg should be 0
    assert_v_tg = abs(vel_aver_opt-v_tgt);    

    %% Reconstruct full gait cycle: starting with right heel strike
    % We reconstruct the full gait cycle from the simulated half gait cycle
    % Identify heel strike
    threshold = 30;
    if exist('HS1','var')
        clear HS1
    end
    % Right heel strike first    
    phase_tran_tgridi = find(GRF_opt_unsc(:,2)<threshold,1,'last');
    if ~isempty(phase_tran_tgridi)        
        if phase_tran_tgridi == N
            temp_idx = find(GRF_opt_unsc(:,2)>threshold,1,'first');
            if ~isempty(temp_idx)
                if temp_idx-1 ~= 0 && ...
                        find(GRF_opt_unsc(temp_idx-1,2)<threshold)
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
    if isempty(phase_tran_tgridi)
        continue;
    end
    % Joint kinematics
    % Qs    
    q_opt_GC                = zeros(N,size(q_opt_unsc.deg,2));    
    q_opt_GC(1:N-IC1i+1,:)  = q_opt_unsc.deg(IC1i:end,:);
    q_opt_GC(N-IC1i+2:N,:)  = q_opt_unsc.deg(1:IC1i-1,:);        
    q_opt_GC(N-IC1i+2:N,jointi.pelvis.tx) = ...
        q_opt_unsc.deg(1:IC1i-1,jointi.pelvis.tx) + ...
        q_opt_unsc_all.deg(end,jointi.pelvis.tx);   
    temp_q_opt_GC_pelvis_tx = q_opt_GC(1,jointi.pelvis.tx);
    q_opt_GC(:,jointi.pelvis.tx) = q_opt_GC(:,jointi.pelvis.tx)-...
        temp_q_opt_GC_pelvis_tx;   

    % Qdots
    qdot_opt_GC = zeros(N,size(q_opt,2));    
    qdot_opt_GC(1:N-IC1i+1,:) = qdot_opt_unsc.deg(IC1i:end,:);
    qdot_opt_GC(N-IC1i+2:N,:) = qdot_opt_unsc.deg(1:IC1i-1,:);

    % Qdotdots
    qdotdot_opt_GC = zeros(N,size(q_opt,2));
    qdotdot_opt_GC(1:N-IC1i+1,:) = qdotdot_opt_unsc.deg(IC1i:end,:);
    qdotdot_opt_GC(N-IC1i+2:N,:) = qdotdot_opt_unsc.deg(1:IC1i-1,:);  

    % Ground reaction forces
    GRF_opt_GC = zeros(N,NGRF);
    GRF_opt_GC(1:N-IC1i+1,:) = GRF_opt_unsc(IC1i:end,:);
    GRF_opt_GC(N-IC1i+2:N,:) = GRF_opt_unsc(1:IC1i-1,:);
    GRF_opt_GC = GRF_opt_GC./(body_weight/100);

    % Joint torques
    tau_opt_GC = zeros(N,size(q_opt,2));    
    tau_opt_GC(1:N-IC1i+1,1:nq.all) = out_res_opt(IC1i:end,1:nq.all);    
    tau_opt_GC(N-IC1i+2:N,:) = out_res_opt(1:IC1i-1,1:nq.all);
    tau_opt_GC = tau_opt_GC./body_mass;    

    % Muscle-Tendon forces
    FTtilde_opt_GC = zeros(N,NMuscle_act);
    FTtilde_opt_GC(1:N-IC1i+1,:) = FTtilde_opt_unsc(IC1i:end,:);
    FTtilde_opt_GC(N-IC1i+2:N,:) = FTtilde_opt_unsc(1:IC1i-1,:);

    % Muscle activations
    a_opt_GC = zeros(N,NMuscle_act);
    a_opt_GC(1:N-IC1i+1,:) = a_opt_unsc(IC1i:end,:);
    a_opt_GC(N-IC1i+2:N,:) = a_opt_unsc(1:IC1i-1,:);

    % Time derivative of muscle-tendon force
    dFTtilde_opt_GC = zeros(N,NMuscle_act);
    dFTtilde_opt_GC(1:N-IC1i+1,:) = dFTtilde_opt_unsc(IC1i:end,:);
    dFTtilde_opt_GC(N-IC1i+2:N,:) = dFTtilde_opt_unsc(1:IC1i-1,:);

    % Muscle excitations
    vA_opt_GC = zeros(N,NMuscle_act);
    vA_opt_GC(1:N-IC1i+1,:) = vA_opt_unsc(IC1i:end,:);
    vA_opt_GC(N-IC1i+2:N,:) = vA_opt_unsc(1:IC1i-1,:);    
    tact = 0.015;
    tdeact = 0.06;
    e_opt_GC = computeExcitationRaasch(a_opt_GC,vA_opt_GC,...
        ones(1,NMuscle_act)*tdeact,ones(1,NMuscle_act)*tact);   
    
    % Arm activations    
    a_a_opt_GC = zeros(N,nq.arms);
    a_a_opt_GC(1:N-IC1i+1,:) = a_a_opt_unsc(IC1i:end,:);
    a_a_opt_GC(N-IC1i+2:N,:) = a_a_opt_unsc(1:IC1i-1,:);
    
    % Arm excitations
    e_a_opt_GC = zeros(N,nq.arms);
    e_a_opt_GC(1:N-IC1i+1,:) = e_a_opt_unsc(IC1i:end,:);
    e_a_opt_GC(N-IC1i+2:N,:) = e_a_opt_unsc(1:IC1i-1,:);

    % Passive joint torques
    Tau_pass_opt_GC = zeros(N,13);
    Tau_pass_opt_GC(1:N-IC1i+1,:) = Tau_pass_opt_all(IC1i:end,:); 
    Tau_pass_opt_GC(N-IC1i+2:N,:) = Tau_pass_opt_all(1:IC1i-1,:); 
    
    % Create .mot file for OpenSim GUI
    q_opt_GUI_GC = zeros(N,1+nq.all+2);
    q_opt_GUI_GC(1:N-IC1i+1,1) = tgrid(:,IC1i:end-1)';    
    q_opt_GUI_GC(N-IC1i+2:N,1)  = tgrid(:,1:IC1i-1)' + tgrid(end);    
    q_opt_GUI_GC(:,2:end-2) = q_opt_GC;
    q_opt_GUI_GC(:,end-1:end) = 1.51*180/pi*ones(N,2); % pro_sup (locked)
    q_opt_GUI_GC(:,1) = q_opt_GUI_GC(:,1)-q_opt_GUI_GC(1,1);
    pathOpenSim = [pathRepo,'/OpenSim'];
    addpath(genpath(pathOpenSim));    
    if writeIKmotion  
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
        % Two gait cycles
        q_opt_GUI_GC_2 = [q_opt_GUI_GC;q_opt_GUI_GC];
        q_opt_GUI_GC_2(N+1:2*N,1) = q_opt_GUI_GC_2(N+1:2*N,1) + ...
            q_opt_GUI_GC_2(end,1) + ...
            q_opt_GUI_GC_2(end,1)-q_opt_GUI_GC_2(end-1,1);
        q_opt_GUI_GC_2(N+1:2*N,jointi.pelvis.tx+1) = ...
            q_opt_GUI_GC_2(N+1:2*N,jointi.pelvis.tx+1) + ...
            q_opt_unsc_all.deg(end,jointi.pelvis.tx);    
        % Muscle activations (to have muscles turning red when activated).
        Acts_opt_GUI = [a_opt_GC;a_opt_GC];
        % Combine data joint angles and muscle activations
        JointAngleMuscleAct.data = [q_opt_GUI_GC_2,Acts_opt_GUI];
        % Get muscle labels
        muscleNamesAll = cell(1,NMuscle-NMuscle_pro);
        for i = 1:NMuscle/2
             muscleNamesAll{i} = [muscleNames{i}(1:end-2),'_l'];
        end 
        for i = 1:NMuscle/2-NMuscle_pro
            muscleNamesAll{i+NMuscle/2} = muscleNames_act{i};
        end        
        % Combine labels joint angles and muscle activations
        JointAngleMuscleAct.labels = JointAngle.labels;
        for i = 1:NMuscle-NMuscle_pro
            JointAngleMuscleAct.labels{i+size(q_opt_GUI_GC_2,2)} = ...
                [muscleNamesAll{i},'/activation'];
        end
        filenameJointAngles = [pathRepo,'/Results/',namescript,...
                '/IK',savename,'.mot'];
        write_motionFile(JointAngleMuscleAct, filenameJointAngles)
    end
    
    %% Reconstruct full gait cycle: starting with left heel strike
    if exist('HS1_l','var')
        clear HS1_l
    end
    % Left heel strike first    
    phase_tran_tgridi_l = find(GRF_opt_unsc(:,5)<threshold,1,'last');
    if ~isempty(phase_tran_tgridi_l)        
        if phase_tran_tgridi_l == N
            temp_idx = find(GRF_opt_unsc(:,5)>threshold,1,'first');
            if ~isempty(temp_idx)
                if temp_idx-1 ~= 0 && ...
                        find(GRF_opt_unsc(temp_idx-1,5)<threshold)
                    phase_tran_tgridi_t = temp_idx;             
                    IC1i_l = phase_tran_tgridi_t;
                    HS1_l = 'l';
                end 
            else            
                IC1i_l = phase_tran_tgridi_l + 1; 
                HS1_l = 'l';
            end
        else            
            IC1i_l = phase_tran_tgridi_l + 1; 
            HS1_l = 'l';
        end        
    end
    if isempty(phase_tran_tgridi_l)
        continue;
    end
    % Joint kinematics
    % Qs     
    q_opt_GC_l                = zeros(N,size(q_opt_unsc.deg,2));    
    q_opt_GC_l(1:N-IC1i_l+1,:)  = q_opt_unsc.deg(IC1i_l:end,:);
    q_opt_GC_l(N-IC1i_l+2:N,:)  = q_opt_unsc.deg(1:IC1i_l-1,:);        
    q_opt_GC_l(N-IC1i_l+2:N,jointi.pelvis.tx) = ...
        q_opt_unsc.deg(1:IC1i_l-1,jointi.pelvis.tx) + ...
        q_opt_unsc_all.deg(end,jointi.pelvis.tx);    
    temp_q_opt_GC_l_pelvis_tx = q_opt_GC_l(1,jointi.pelvis.tx);
    q_opt_GC_l(:,jointi.pelvis.tx) = q_opt_GC_l(:,jointi.pelvis.tx)-...
        temp_q_opt_GC_l_pelvis_tx;    

    % Qdots
    qdot_opt_GC_l = zeros(N,size(q_opt,2));    
    qdot_opt_GC_l(1:N-IC1i_l+1,:) = qdot_opt_unsc.deg(IC1i_l:end,:);
    qdot_opt_GC_l(N-IC1i_l+2:N,:) = qdot_opt_unsc.deg(1:IC1i_l-1,:);

    % Qdotdots
    qdotdot_opt_GC_l = zeros(N,size(q_opt,2));
    qdotdot_opt_GC_l(1:N-IC1i_l+1,:) = qdotdot_opt_unsc.deg(IC1i_l:end,:);
    qdotdot_opt_GC_l(N-IC1i_l+2:N,:) = qdotdot_opt_unsc.deg(1:IC1i_l-1,:);  

    % Ground reaction forces
    GRF_opt_GC_l = zeros(N,NGRF);
    GRF_opt_GC_l(1:N-IC1i_l+1,:) = GRF_opt_unsc(IC1i_l:end,:);
    GRF_opt_GC_l(N-IC1i_l+2:N,:) = GRF_opt_unsc(1:IC1i_l-1,:);
    GRF_opt_GC_l = GRF_opt_GC_l./(body_weight/100);
    
    % Joint torques
    tau_opt_GC_l = zeros(N,size(q_opt,2));    
    tau_opt_GC_l(1:N-IC1i_l+1,1:nq.all) = out_res_opt(IC1i_l:end,1:nq.all);    
    tau_opt_GC_l(N-IC1i_l+2:N,:) = out_res_opt(1:IC1i_l-1,1:nq.all);
    tau_opt_GC_l = tau_opt_GC_l./body_mass;    

    % Muscle-Tendon forces
    FTtilde_opt_GC_l = zeros(N,NMuscle_act);
    FTtilde_opt_GC_l(1:N-IC1i_l+1,:) = FTtilde_opt_unsc(IC1i_l:end,:);
    FTtilde_opt_GC_l(N-IC1i_l+2:N,:) = FTtilde_opt_unsc(1:IC1i_l-1,:);

    % Muscle activations
    a_opt_GC_l = zeros(N,NMuscle_act);
    a_opt_GC_l(1:N-IC1i_l+1,:) = a_opt_unsc(IC1i_l:end,:);
    a_opt_GC_l(N-IC1i_l+2:N,:) = a_opt_unsc(1:IC1i_l-1,:);

    % Time derivative of muscle-tendon force
    dFTtilde_opt_GC_l = zeros(N,NMuscle_act);
    dFTtilde_opt_GC_l(1:N-IC1i_l+1,:) = dFTtilde_opt_unsc(IC1i_l:end,:);
    dFTtilde_opt_GC_l(N-IC1i_l+2:N,:) = dFTtilde_opt_unsc(1:IC1i_l-1,:);

    % Muscle excitations
    vA_opt_GC_l = zeros(N,NMuscle_act);
    vA_opt_GC_l(1:N-IC1i_l+1,:) = vA_opt_unsc(IC1i_l:end,:);
    vA_opt_GC_l(N-IC1i_l+2:N,:) = vA_opt_unsc(1:IC1i_l-1,:);
    tact = 0.015;
    tdeact = 0.06;
    e_opt_GC_l = computeExcitationRaasch(a_opt_GC_l,vA_opt_GC_l,...
        ones(1,NMuscle_act)*tdeact,ones(1,NMuscle_act)*tact);   
    
    % Arm activations    
    a_a_opt_GC_l = zeros(N,nq.arms);
    a_a_opt_GC_l(1:N-IC1i_l+1,:) = a_a_opt_unsc(IC1i_l:end,:);
    a_a_opt_GC_l(N-IC1i_l+2:N,:) = a_a_opt_unsc(1:IC1i_l-1,:);
    
    % Arm excitations
    e_a_opt_GC_l = zeros(N,nq.arms);
    e_a_opt_GC_l(1:N-IC1i_l+1,:) = e_a_opt_unsc(IC1i_l:end,:);
    e_a_opt_GC_l(N-IC1i_l+2:N,:) = e_a_opt_unsc(1:IC1i_l-1,:);

    % Passive joint torques
    Tau_pass_opt_GC_l = zeros(N,13);
    Tau_pass_opt_GC_l(1:N-IC1i_l+1,:) = Tau_pass_opt_all(IC1i_l:end,:); 
    Tau_pass_opt_GC_l(N-IC1i_l+2:N,:) = Tau_pass_opt_all(1:IC1i_l-1,:); 
    
    % Create .mot file for OpenSim GUI
    q_opt_GUI_GC_l = zeros(N,1+nq.all+2);
    q_opt_GUI_GC_l(1:N-IC1i_l+1,1) = tgrid(:,IC1i_l:end-1)';    
    q_opt_GUI_GC_l(N-IC1i_l+2:N,1)  = tgrid(:,1:IC1i_l-1)' + tgrid(end);    
    q_opt_GUI_GC_l(:,2:end-2) = q_opt_GC_l;
    q_opt_GUI_GC_l(:,end-1:end) = 1.51*180/pi*ones(N,2); % pro_sup (locked)
    q_opt_GUI_GC_l(:,1) = q_opt_GUI_GC_l(:,1)-q_opt_GUI_GC_l(1,1);
    pathOpenSim = [pathRepo,'/OpenSim'];
    addpath(genpath(pathOpenSim));       
    if writeIKmotion  
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
        % Two gait cycles
        q_opt_GUI_GC_l_2 = [q_opt_GUI_GC_l;q_opt_GUI_GC_l];
        q_opt_GUI_GC_l_2(N+1:2*N,1) = q_opt_GUI_GC_l_2(N+1:2*N,1) + ...
            q_opt_GUI_GC_l_2(end,1) + ...
            q_opt_GUI_GC_l_2(end,1)-q_opt_GUI_GC_l_2(end-1,1);
        q_opt_GUI_GC_l_2(N+1:2*N,jointi.pelvis.tx+1) = ...
            q_opt_GUI_GC_l_2(N+1:2*N,jointi.pelvis.tx+1) + ...
            q_opt_unsc_all.deg(end,jointi.pelvis.tx);   
        % Muscle activations (to have muscles turning red when activated).
        Acts_opt_l_GUI = [a_opt_GC_l;a_opt_GC_l];
        % Combine data joint angles and muscle activations
        JointAngleMuscleAct_l.data = [q_opt_GUI_GC_l_2,Acts_opt_l_GUI];
        % Get muscle labels
        muscleNamesAll = cell(1,NMuscle-NMuscle_pro);
        for i = 1:NMuscle/2
             muscleNamesAll{i} = [muscleNames{i}(1:end-2),'_l'];
        end 
        for i = 1:NMuscle/2-NMuscle_pro
            muscleNamesAll{i+NMuscle/2} = muscleNames_act{i};
        end        
        % Combine labels joint angles and muscle activations
        JointAngleMuscleAct_l.labels = JointAngle.labels;
        for i = 1:NMuscle-NMuscle_pro
            JointAngleMuscleAct_l.labels{i+size(q_opt_GUI_GC_2,2)} = ...
                [muscleNamesAll{i},'/activation'];
        end
        filenameJointAngles = [pathRepo,'/Results/',namescript,...
                '/IK',savename,'_l.mot'];
        write_motionFile(JointAngleMuscleAct_l, filenameJointAngles)
    end
        
    %% Metabolic cost of transport
    q_opt_GC_rad = q_opt_GC;
    q_opt_GC_rad(:,[1:3,7:end]) = q_opt_GC_rad(:,[1:3,7:end]).*pi/180;
    qdot_opt_GC_rad = qdot_opt_GC;
    qdot_opt_GC_rad(:,[1:3,7:end])= qdot_opt_GC_rad(:,[1:3,7:end]).*pi/180;    
    % Pre-allocations
    e_mo_opt        = zeros(N, 1);
    vMtilde_opt_all = zeros(N, NMuscle_act);
    lMtilde_opt_all = zeros(N, NMuscle_act);
    for nn = 1:N
        % Get muscle-tendon lengths, velocities, moment arms
        % Left leg
        qin_l_opt = [q_opt_GC_rad(nn,jointi.hip_flex.l),...
            q_opt_GC_rad(nn,jointi.hip_add.l), ...
            q_opt_GC_rad(nn,jointi.hip_rot.l), ...
            q_opt_GC_rad(nn,jointi.knee.l), ...
            q_opt_GC_rad(nn,jointi.ankle.l),...
            q_opt_GC_rad(nn,jointi.subt.l),...
            q_opt_GC_rad(nn,jointi.trunk.ext),...
            q_opt_GC_rad(nn,jointi.trunk.ben),...
            q_opt_GC_rad(nn,jointi.trunk.rot)];  
        qdotin_l_opt = [qdot_opt_GC_rad(nn,jointi.hip_flex.l),...
            qdot_opt_GC_rad(nn,jointi.hip_add.l),...
            qdot_opt_GC_rad(nn,jointi.hip_rot.l),...
            qdot_opt_GC_rad(nn,jointi.knee.l),...
            qdot_opt_GC_rad(nn,jointi.ankle.l),...
            qdot_opt_GC_rad(nn,jointi.subt.l),...
            qdot_opt_GC_rad(nn,jointi.trunk.ext),...
            qdot_opt_GC_rad(nn,jointi.trunk.ben),...
            qdot_opt_GC_rad(nn,jointi.trunk.rot)];  
        [lMTk_l_opt,vMTk_l_opt,~] = ...
            f_lMT_vMT_dM(qin_l_opt,qdotin_l_opt);    
        % Right leg
        qin_r_opt = [q_opt_GC_rad(nn,jointi.hip_flex.r),...
            q_opt_GC_rad(nn,jointi.hip_add.r),...
            q_opt_GC_rad(nn,jointi.hip_rot.r),...
            q_opt_GC_rad(nn,jointi.knee.r),...
            q_opt_GC_rad(nn,jointi.ankle.r),...
            q_opt_GC_rad(nn,jointi.subt.r),...
            q_opt_GC_rad(nn,jointi.trunk.ext),...
            q_opt_GC_rad(nn,jointi.trunk.ben),...
            q_opt_GC_rad(nn,jointi.trunk.rot)];  
        qdotin_r_opt = [qdot_opt_GC_rad(nn,jointi.hip_flex.r),...
            qdot_opt_GC_rad(nn,jointi.hip_add.r),...
            qdot_opt_GC_rad(nn,jointi.hip_rot.r),...
            qdot_opt_GC_rad(nn,jointi.knee.r),...
            qdot_opt_GC_rad(nn,jointi.ankle.r),...
            qdot_opt_GC_rad(nn,jointi.subt.r),...
            qdot_opt_GC_rad(nn,jointi.trunk.ext),...
            qdot_opt_GC_rad(nn,jointi.trunk.ben),...
            qdot_opt_GC_rad(nn,jointi.trunk.rot)];      
        [lMTk_r_opt,vMTk_r_opt,~] = ...
            f_lMT_vMT_dM(qin_r_opt,qdotin_r_opt);
        % Both legs
        lMTk_lr_opt = [lMTk_l_opt([1:43,47:49],1);lMTk_r_opt(musi_pro,1)];
        vMTk_lr_opt = [vMTk_l_opt([1:43,47:49],1);vMTk_r_opt(musi_pro,1)];   
        [Hilldiff_optt,FT_optt,Fce_optt,Fpass_optt,Fiso_optt,...
            vMmax_optt,massM_optt] = ...
                f_forceEquilibrium_FtildeState_all_tendon(...
                a_opt_GC(nn,:)',FTtilde_opt_GC(nn,:)',...
                dFTtilde_opt_GC(nn,:)',full(lMTk_lr_opt),...
                full(vMTk_lr_opt),tensions,aTendon,shift);  
        [~,lMtilde_opt] = f_FiberLength_TendonForce_tendon(...
            FTtilde_opt_GC(nn,:)',full(lMTk_lr_opt),aTendon,shift);         
        lMtilde_opt_all(nn,:) = full(lMtilde_opt)';        
        [vM_opt,vMtilde_opt] = ...
            f_FiberVelocity_TendonForce_tendon(FTtilde_opt_GC(nn,:)',...
            dFTtilde_opt_GC(nn,:)',full(lMTk_lr_opt),full(vMTk_lr_opt),...
            aTendon,shift);        
        vMtilde_opt_all(nn,:) = full(vMtilde_opt)';  
        % Bhargava et al. (2004)
        [~,~,~,~,~,e_mot] = ...
            fgetMetabolicEnergySmooth2004all(a_opt_GC(nn,:)',...
            a_opt_GC(nn,:)',full(lMtilde_opt),full(vM_opt),...
            full(Fce_optt),full(Fpass_optt),full(massM_optt),pctsts,...
            full(Fiso_optt),MTparameters_m(1,:)',body_mass,10);
        e_mo_opt(nn,:) = full(e_mot)';         
    end               
    % Get COT
    dist_trav_opt_GC = q_opt_GC_rad(end,jointi.pelvis.tx) - ...
        q_opt_GC_rad(1,jointi.pelvis.tx);  
    time_GC = q_opt_GUI_GC(:,1);
    e_mo_opt_tr = trapz(time_GC,e_mo_opt);
    % Cost of transport: J/kg/m
    % Energy model from Bhargava et al. (2004)
    COT_opt = e_mo_opt_tr/body_mass/dist_trav_opt_GC;   
    
    %% Optimal cost and CPU time
    pathDiary = [pathresults,'/',namescript,'/D',savename];
    [CPU_IPOPT,CPU_NLP,~,Cost,~,~,~,~,OptSol] = readDiary(pathDiary);
    
    %% Save results       
    if saveResults
        if (exist([pathresults,'/',namescript,...
                '/Results_prosthesis.mat'],'file')==2) 
            load([pathresults,'/',namescript,'/Results_prosthesis.mat']);
        else
            Results_prosthesis.(['Speed_',num2str(v_tgt_id*100)]). ...  
                (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
                (['W_MuscleActivity_',num2str(W.A)]). ...            
                (['W_JointAcceleration_',num2str(W.Ak)]). ...
                (['W_PassiveTorque_',num2str(W.passMom)]). ...
                (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
                (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
                (['InitialGuessType_',num2str(IGsel)]). ...
                (['InitialGuessMode_',num2str(IGm)]). ...
                (['InitialGuessCase_',num2str(IGcase)]). ... 
                (['WeaknessHipActuators_',num2str(h_weak)]). ...
                (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
                (['MetabolicEnergyModel_',num2str(mE)]). ...            
                (['ContactModel_',num2str(cm)]). ...  
                (['Number_MeshIntervals_',num2str(N)]). ...
                (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
                (['Stiffness_',num2str(kstiff)]) = struct('Qs_opt_r',[]);
        end    
        % Put data in structure
        Results_prosthesis.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['Stiffness_',num2str(kstiff)]).Qs_opt_r = q_opt_GC;
        Results_prosthesis.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['Stiffness_',num2str(kstiff)]).Qs_opt_l = q_opt_GC_l;
        Results_prosthesis.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['Stiffness_',num2str(kstiff)]).Qdots_opt_r = qdot_opt_GC;
        Results_prosthesis.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['Stiffness_',num2str(kstiff)]).Qdots_opt_l = qdot_opt_GC_l;
        Results_prosthesis.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['Stiffness_',num2str(kstiff)]).GRFs_opt_r = GRF_opt_GC;
        Results_prosthesis.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['Stiffness_',num2str(kstiff)]).GRFs_opt_l = GRF_opt_GC_l;
        Results_prosthesis.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['Stiffness_',num2str(kstiff)]).Ts_opt_r = tau_opt_GC;
        Results_prosthesis.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['Stiffness_',num2str(kstiff)]).Ts_opt_l = tau_opt_GC_l;
        Results_prosthesis.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['Stiffness_',num2str(kstiff)]).Acts_opt = a_opt_GC;
        Results_prosthesis.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['Stiffness_',num2str(kstiff)]).Acts_opt_l = a_opt_GC_l;        
        Results_prosthesis.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['Stiffness_',num2str(kstiff)]).COT_opt = COT_opt;
        Results_prosthesis.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['Stiffness_',num2str(kstiff)]).CPU_IPOPT = CPU_IPOPT;
        Results_prosthesis.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['Stiffness_',num2str(kstiff)]).CPU_NLP = CPU_NLP;
        Results_prosthesis.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['Stiffness_',num2str(kstiff)]).Cost = Cost;
        Results_prosthesis.(['Speed_',num2str(v_tgt_id*100)]). ...  
            (['W_MetabolicEnergyRate_',num2str(W.E)]). ...
            (['W_MuscleActivity_',num2str(W.A)]). ...            
            (['W_JointAcceleration_',num2str(W.Ak)]). ...
            (['W_PassiveTorque_',num2str(W.passMom)]). ...
            (['W_ArmExcitation_',num2str(W.ArmE)]). ...               
            (['Power_MetabolicEnergyRate_',num2str(exp_E)]). ...
            (['InitialGuessType_',num2str(IGsel)]). ...
            (['InitialGuessMode_',num2str(IGm)]). ...
            (['InitialGuessCase_',num2str(IGcase)]). ... 
            (['WeaknessHipActuators_',num2str(h_weak)]). ...
            (['WeaknessAnklePlantarflexors_',num2str(pf_weak)]). ...   
            (['MetabolicEnergyModel_',num2str(mE)]). ...            
            (['ContactModel_',num2str(cm)]). ...  
            (['Number_MeshIntervals_',num2str(N)]). ...
            (['MaximumContractionVelocity_',num2str(vMax_s)]). ...
            (['Stiffness_',num2str(kstiff)]).OptSol = OptSol;
        Results_prosthesis.colheaders.joints = joints;
        Results_prosthesis.colheaders.GRF = {'fore_aft_r','vertical_r',...
        'lateral_r','fore_aft_l','vertical_l','lateral_l'};
        for i = 1:NMuscle/2
            Results_prosthesis.colheaders.muscles{i} = ...
                [muscleNames{i}(1:end-2),'_l'];
        end 
        for i = 1:NMuscle/2-NMuscle_pro
            Results_prosthesis.colheaders.muscles{i+NMuscle/2} = ...
                muscleNames_act{i};
        end        
        % Save data
        save([pathresults,'/',namescript,'/Results_prosthesis.mat'],...
            'Results_prosthesis');
    end
end
end