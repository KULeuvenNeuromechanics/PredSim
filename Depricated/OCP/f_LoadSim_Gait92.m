function [R] = f_LoadSim_Gait92(ResultsFolder,loadname,varargin)


%% Notes

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

if strcmp(loadname(end-3:end),'.mat')
    loadname = loadname(1:end-4);
end

pathmain = mfilename('fullpath');
[filepath,~,~] =  fileparts(pathmain);
[pathRepo,~,~] = fileparts(filepath);
OutFolder = fullfile(pathRepo,'Results',ResultsFolder);
Outname = fullfile(OutFolder,[loadname '.mat']);
load(Outname,'w_opt','stats','Sopt','ExoVect','setup');
S = Sopt;

body_mass = S.mass;
%% User inputs (typical settings structure)
% load default CasadiFunctions

% flow control
writeIKmotion   = 1; % set to 1 to write .mot file

% settings for optimization
v_tgt       = S.v_tgt;      % average speed
N           = S.N;          % number of mesh intervals
W           = S.W;
exp_E       = S.W.exp_E;    % power metabolic energy

% ipopt options
tol_ipopt       = S.tol_ipopt;

%% Settings

import casadi.*
if ~isfield(S,'subject') || isempty(S.subject)
    S.subject = 'subject1';
end
subject = S.subject;

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
roti                = [jointi.pelvis.tilt:jointi.pelvis.rot,...
    jointi.hip_flex.l:jointi.elb.r];
% Number of degrees of freedom for later use
nq.all      = length(residualsi); % all
nq.abs      = length(ground_pelvisi); % ground-pelvis
nq.trunk    = length(trunki); % trunk
nq.arms     = length(armsi); % arms
nq.mtp     = length(mtpi); % arms
nq.leg      = 10; % #joints needed for polynomials

%% Load external functions
% The external function performs inverse dynamics through the
% OpenSim/Simbody C++ API. This external function is compiled as a dll from
% which we create a Function instance using CasADi in MATLAB. More details
% about the external function can be found in the documentation.
pathmain = pwd;
% We use different external functions, since we also want to access some
% parameters of the model in a post-processing phase.
pathExternalFunctions = [pathRepo,'/ExternalFunctions'];
% Loading external functions.
cd(pathExternalFunctions);

% external function used in the optimization
ExtF = S.ExternalFunc;
F = external('F',ExtF);

% external function selected for post-processing
if ~isempty(varargin)
    S.ExternalFunc2 = varargin{1};
end

if strcmp(S.ExternalFunc2,'Schertzer2014_pp.ddl')
    S.ExternalFunc2 = 'Schertzer2014_pp.dll';
end
if strcmp(S.ExternalFunc2,'Analyse_s1Pog.ddl')
    S.ExternalFunc2 = 'Analyse_s1Pog.dll';
end

if isfield(S,'ExternalFunc2')    
    F1 = external('F',S.ExternalFunc2);
elseif F.nnz_in ==  nq.all*3
    warning(['No external was selected for post processing. Post processing',...
        'is done with Browning_2008_pp.dll']);
    F1 = external('F', 'Browning_2008_pp.dll');
elseif F.nnz_in ==  nq.all*3 + 2
    F1 = external('F','SimExo_3D_ExportAll.dll');
    disp('Selected SimExo_3D_ExportAll.dll as external function because S.ExternalFunc2 was not specified');
end
cd(pathmain);



%% Test the type of simulation
ExoImplementation = 'IdealAnkle'; % (1) ankle actuation (2) passive afo, (3) torque actuator
if F.nnz_in == nq.all*3+2
    % After the first simulations report (May 12 2020), we changed the implemenation of the exoskeleton assistance
    % to a torque actuator at the calcaneus and tibia. We also changed from two .dll files (one for optimization and one
    % for post processing) to one .dll file for everything.
    ExoImplementation = 'TorqueTibiaCalcn';
elseif strcmp(ExtF,'PredSim_3D_GRF.dll') && isfield(S,'AFO_stiffness')
    % AFO modelled as in Nuckols 2019: Ultrasound imaging links soleus
    % muscle neuromechanics and energetics during human walking with
    % elastic ankle exoskeletons
    disp('Post processing of passive AFO as in Nuckols - ultrasound paper');
    ExoImplementation = 'Nuckols2019';
elseif F.nnz_in == nq.all*3+6
    % typtical implementation when simulation with ankle-knee-hip
    % exoskeleton
    ExoImplementation = 'TorqueTibiaCalcn';
else
    disp('Default processing with active exoskeleton (Poggensee paper)')
end




%% Second, origins bodies.
% Calcaneus
if strcmp(ExoImplementation,'IdealAnkle') || strcmp(ExoImplementation,'Nuckols2019')
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
end
% adapt indexes GRF in Nuckols simulations
if strcmp(ExoImplementation,'Nuckols2019')
    GRFi.r      = GRFi.r + 20;
    GRFi.l      = GRFi.l + 20;
    GRFi.all    = [GRFi.r,GRFi.l];
end

if strcmp(ExoImplementation,'TorqueTibiaCalcn') || F1.nnz_out == 73
    calcOr.r    = [38 40];
    calcOr.l    = [41 43];
    calcOr.all  = [calcOr.r,calcOr.l];
    NcalcOr     = length(calcOr.all);
    % Femurs
    femurOr.r   = [44 46];
    femurOr.l   = [47 49];
    femurOr.all = [femurOr.r,femurOr.l];
    NfemurOr    = length(femurOr.all);
    % Hands
    handOr.r    = [50 52];
    handOr.l    = [53 55];
    handOr.all  = [handOr.r,handOr.l];
    NhandOr     = length(handOr.all);
    % Tibias
    tibiaOr.r   = [56 58];
    tibiaOr.l   = [59 61];
    tibiaOr.all = [tibiaOr.r,tibiaOr.l];
    NtibiaOr    = length(tibiaOr.all);
    % External function: F1 (post-processing purpose only)
    % Ground reaction forces (GRFs)
    GRFi.r      = 32:34;
    GRFi.l      = 35:37;
    GRFi.all    = [GRFi.r,GRFi.l];
    NGRF        = length(GRFi.all);
    % toe joints
    toesOr.r   = [62 64];
    toesOr.l   = [65 67];
    toesOr.all = [toesOr.r,toesOr.l];
    NtoesOr    = length(toesOr.all);
    
    % Origins calcaneus (3D)
    calcOrall.r     = 38:40;
    calcOrall.l     = 41:43;
    calcOrall.all   = [calcOrall.r,calcOrall.l];
    NcalcOrall      = length(calcOrall.all);
    
    % COP information
    GroundT.r = [68 69 70];
    GroundT.l = [71 72 73];
end


%% CasADi functions
% We create several CasADi functions for later use
pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
PathDefaultFunc = fullfile(pathCasADiFunctions,S.CasadiFunc_Folders);
cd(PathDefaultFunc);
f_FiberLength_TendonForce_tendon = Function.load('f_FiberLength_TendonForce_tendon');
f_FiberVelocity_TendonForce_tendon = Function.load('f_FiberVelocity_TendonForce_tendon');
f_forceEquilibrium_FtildeState_all_tendon = Function.load('f_forceEquilibrium_FtildeState_all_tendon');
f_J2    = Function.load('f_J2');
f_J23   = Function.load('f_J23');
f_J25   = Function.load('f_J25');
f_J8    = Function.load('f_J8');
f_J92   = Function.load('f_J92');
f_J92exp = Function.load('f_J92exp');
f_Jnn3  = Function.load('f_Jnn3');
f_lMT_vMT_dM = Function.load('f_lMT_vMT_dM');
f_AllPassiveTorques = Function.load('f_AllPassiveTorques');
fgetMetabolicEnergySmooth2004all = Function.load('fgetMetabolicEnergySmooth2004all');
cd(pathmain);

%% load the metalbolic energy equations
PathEnergyEq = fullfile(pathCasADiFunctions,'EnergyModels');
cd(PathEnergyEq);
fgetMetabolicEnergySmooth2003all    = Function.load('fgetMetabolicEnergySmooth2003all');
fgetMetabolicEnergySmooth2010all    = Function.load('fgetMetabolicEnergySmooth2010all');
fgetMetabolicEnergySmooth2016all    = Function.load('fgetMetabolicEnergySmooth2016all');
fgetMetabolicEnergySmooth2010all_hl = Function.load('fgetMetabolicEnergySmooth2010all_hl');
fgetMetabolicEnergySmooth2010all_neg= Function.load('fgetMetabolicEnergySmooth2010all_neg');
fgetMetabolicEnergy_MargariaSmooth  = Function.load('fgetMetabolicEnergy_MargariaSmooth');
cd(pathmain);

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

% indices of muscles
musi = MuscleIndices(muscleNames(1:end-3),muscleNames);
NMuscle = length(muscleNames(1:end-3))*2;

% Muscle indices for later use
File_MTparameters = fullfile(PathDefaultFunc,'MTparameters.mat');
if exist(File_MTparameters,'file')
    load(File_MTparameters,'MTparameters');
else
    % This file was saved in another location in the old implementation
    pathmusclemodel = fullfile(pathRepo,'MuscleModel',subject);
    ExtPoly = '_mtp';
    load([pathmusclemodel,'/MTparameters_',subject, ExtPoly, '.mat'],'MTparameters');
end
MTparameters_m = [MTparameters(:,musi),MTparameters(:,musi)];

% path to the polynomial functions
if isfield(S,'PolyFolder') && ~isempty(S.PolyFolder)
    % default location
    pathpolynomial = fullfile(pathRepo,'Polynomials',S.PolyFolder);
else
    % old version (we still want to be able to process these results)
    pathpolynomial = fullfile(pathRepo,'Polynomials',S.subject);
end
PolyFile = [pathpolynomial,'/muscle_spanning_joint_INFO_',subject,'_mtp.mat'];
if exist(PolyFile,'file')
    tl = load(PolyFile);
else
    PolyFile = [pathpolynomial,'/muscle_spanning_joint_INFO.mat'];
    tl = load(PolyFile);
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


%% file with mass of muscles
MassFile = fullfile(PathDefaultFunc,'MassM.mat');
if exist(MassFile,'file')
    MuscleMass = load(MassFile);
else
    MassFile = fullfile(pathCasADiFunctions,'MassM.mat');
    MuscleMass =load(MassFile);
end

%% Exoskeleton torque (needed to process old simulation files)
if ~exist('ExoVect','var')
    load(Outname,'ExoControl');
    if ~isempty(ExoControl)
        ExoVect = [ExoControl.Tankle_l; ExoControl.Tankle_r];
    else
        ExoVect = zeros(2,N);
    end
end
%% Joints
joints = {'pelvis_tilt','pelvis_list','pelvis_rotation','pelvis_tx',...
    'pelvis_ty','pelvis_tz','hip_flexion_l','hip_adduction_l',...
    'hip_rotation_l','hip_flexion_r','hip_adduction_r','hip_rotation_r',...
    'knee_angle_l','knee_angle_r','ankle_angle_l','ankle_angle_r',...
    'subtalar_angle_l','subtalar_angle_r','mtp_angle_l','mtp_angle_r',...
    'lumbar_extension','lumbar_bending','lumbar_rotation','arm_flex_l',...
    'arm_add_l','arm_rot_l','arm_flex_r','arm_add_r','arm_rot_r',...
    'elbow_flex_l','elbow_flex_r'};

%% get scaling
scaling = setup.scaling;

%% Index helpers

% indexes to select kinematics left and right leg
IndexLeft = [jointi.hip_flex.l jointi.hip_add.l jointi.hip_rot.l, ...
    jointi.knee.l jointi.ankle.l jointi.subt.l jointi.mtp.l,...
    jointi.trunk.ext, jointi.trunk.ben, jointi.trunk.rot];
IndexRight = [jointi.hip_flex.r jointi.hip_add.r jointi.hip_rot.r, ...
    jointi.knee.r jointi.ankle.r jointi.subt.r jointi.mtp.r,...
    jointi.trunk.ext, jointi.trunk.ben, jointi.trunk.rot];

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


%% Read from the vector with optimization results

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
% when optimizing exoskeleton assistance
if isfield(S,'OptTexo_Ankle')
    if S.OptTexo_Ankle.Bool
        ExoVect =reshape(w_opt(starti:starti+N*2-1),2,N);
        starti = starti + 2*N;
    elseif S.OptTexo_AnkleKneeHip.Bool
        ExoVect =reshape(w_opt(starti:starti+N*6-1),6,N);
        starti = starti + 6*N;
    end
end
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
Foutk_opt                   = zeros(N,F1.nnz_out);
Tau_passk_opt_all           = zeros(N,nq.all-nq.abs);
if S.ExoBool == 1 && strcmp(ExoImplementation,'TorqueTibiaCalcn')
    Foutk_opt_Exo         = zeros(N,F1.nnz_out);
end
ExoZeroT = zeros(length(ExoVect(:,1)),1);

for i = 1:N
    % ID moments
    if F.nnz_in > nq.all*3
        if S.ExoBool == 1 && strcmp(ExoImplementation,'TorqueTibiaCalcn')
            % compute torque with exoskeleton support
            [res2] = F1([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)'; -ExoVect(:,i)]);
            [res2_or] = F([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)'; -ExoVect(:,i)]); % ext func used in optimization
        end
        % compute torque without exoskeleton support
        [res] = F1([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)'; ExoZeroT]);
        [res_or] = F([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)';ExoZeroT]); % ext func used in optimization
    else
        [res] = F1([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)']);
        [res_or] = F([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)']); % ext func used in optimization
    end
    Foutk_opt(i,:) = full(res);
    Foutk_opt(i,1:nq.all) = full(res_or(1:nq.all)); % extract ID moments from external function used in the optimization
    if S.ExoBool == 1 && strcmp(ExoImplementation,'TorqueTibiaCalcn')
        Foutk_opt_Exo(i,:) = full(res2);
        Foutk_opt_Exo(i,1:nq.all) = full(res2_or(1:nq.all));% ID moments based on original function (just to be sure)
    end
    % passive moments
    Tau_passk_opt_all(i,:) = full(f_AllPassiveTorques(q_opt_unsc_all.rad(i+1,:),qdot_opt_unsc_all.rad(i+1,:)));
end
GRFk_opt = Foutk_opt(:,GRFi.all);



%% Joint torques and ground reaction forces at collocation points
Xj_Qs_Qdots_opt             = zeros(d*N,2*nq.all);
Xj_Qs_Qdots_opt(:,1:2:end)  = q_col_opt_unsc.rad;
Xj_Qs_Qdots_opt(:,2:2:end)  = qdot_col_opt_unsc.rad;
Xj_Qdotdots_opt             = qdotdot_col_opt_unsc.rad;
Foutj_opt                   = zeros(d*N,F1.nnz_out);
Tau_passj_opt_all           = zeros(d*N,nq.all-nq.abs);
if S.ExoBool == 1 && strcmp(ExoImplementation,'TorqueTibiaCalcn')
    Foutj_opt_Exo         = zeros(d*N,F1.nnz_out);
end
for i = 1:d*N
    iMesh = ceil(i/d);
    % inverse dynamics
    if F1.nnz_in > nq.all*3
        if S.ExoBool == 1 && strcmp(ExoImplementation,'TorqueTibiaCalcn')
            [res2] = F1([Xj_Qs_Qdots_opt(i,:)';Xj_Qdotdots_opt(i,:)'; -ExoVect(:,iMesh)]);
            [res2_or] = F([Xj_Qs_Qdots_opt(i,:)';Xj_Qdotdots_opt(i,:)'; -ExoVect(:,iMesh)]);
        end
        [res] = F1([Xj_Qs_Qdots_opt(i,:)';Xj_Qdotdots_opt(i,:)'; ExoZeroT]);
        [res_or] = F([Xj_Qs_Qdots_opt(i,:)';Xj_Qdotdots_opt(i,:)'; ExoZeroT]);
    else
        [res] = F1([Xj_Qs_Qdots_opt(i,:)';Xj_Qdotdots_opt(i,:)']);
        [res_or] = F([Xj_Qs_Qdots_opt(i,:)';Xj_Qdotdots_opt(i,:)']);
    end
    Foutj_opt(i,:) = full(res);
    Foutj_opt(i,1:nq.all) = full(res_or(1:nq.all)); % extract ID moments from external function used in the optimization
    if S.ExoBool == 1 && strcmp(ExoImplementation,'TorqueTibiaCalcn')
        Foutj_opt_Exo(i,:) = full(res2);
        Foutj_opt_Exo(i,1:nq.all) = full(res2_or(1:nq.all));% ID moments based on original function (just to be sure)
    end
    % passive torques
    Tau_passj_opt_all(i,:) = full(f_AllPassiveTorques(q_col_opt_unsc.rad(i,:),qdot_col_opt_unsc.rad(i,:)));
end
Tau_passj_J = Tau_passj_opt_all(:,[1:12 15:end]);


%% Stride length and width
% For the stride length we also need the values at the end of the
% interval so N+1 where states but not controls are defined
Xk_Qs_Qdots_opt_all = zeros(N+1,2*size(q_opt_unsc_all.rad,2));
Xk_Qs_Qdots_opt_all(:,1:2:end)  = q_opt_unsc_all.rad;
Xk_Qs_Qdots_opt_all(:,2:2:end)  = qdot_opt_unsc_all.rad;
% We just want to extract the positions of the calcaneus origins so we
% do not really care about Qdotdot that we set to 0
Xk_Qdotdots_opt_all = zeros(N+1,size(q_opt_unsc_all.rad,2));
out_res_opt_all = zeros(N+1,F1.nnz_out);
ndof = size(q_opt_unsc_all.rad,2);
for i = 1:N+1
    if F1.nnz_in == nq.all*3
        [res] = F1([Xk_Qs_Qdots_opt_all(i,:)';Xk_Qdotdots_opt_all(i,:)']);
        [res_or] = F([Xk_Qs_Qdots_opt_all(i,:)';Xk_Qdotdots_opt_all(i,:)']);
    else
        [res] = F1([Xk_Qs_Qdots_opt_all(i,:)';Xk_Qdotdots_opt_all(i,:)'; ExoZeroT]);
        [res_or] = F([Xk_Qs_Qdots_opt_all(i,:)';Xk_Qdotdots_opt_all(i,:)'; ExoZeroT]);
    end
    out_res_opt_all(i,:) = full(res);
    out_res_opt_all(i,1:nq.all) = full(res_or(1:nq.all)); % ID moments based on original function (just to be sure)
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
        lMTk_lr_opt_all = [lMTk_l_opt_all([1:43,47:49],1);lMTk_r_opt_all(1:46,1)];
        vMTk_lr_opt_all = [vMTk_l_opt_all([1:43,47:49],1);vMTk_r_opt_all(1:46,1)];
        % force equilibirum
        [~,~,Fce_opt_all,Fpass_opt_all,Fiso_opt_all] = ...
            f_forceEquilibrium_FtildeState_all_tendon(...
            a_col_opt_unsc(count,:)',FTtilde_col_opt_unsc(count,:)',...
            dFTtilde_col_opt_unsc(count,:)',full(lMTk_lr_opt_all),...
            full(vMTk_lr_opt_all),tensions);
        % muscle-tendon kinematics
        [~,lMtilde_opt_all] = f_FiberLength_TendonForce_tendon(...
            FTtilde_col_opt_unsc(count,:)',full(lMTk_lr_opt_all));
        [vM_opt_all,~] = f_FiberVelocity_TendonForce_tendon(...
            FTtilde_col_opt_unsc(count,:)',...
            dFTtilde_col_opt_unsc(count,:)',full(lMTk_lr_opt_all),...
            full(vMTk_lr_opt_all));
        
        % Bhargava et al. (2004)
        if isfield(S,'EModel')
            if strcmp(S.EModel,'Bhargava2014')
                [e_tot_all,~,~,~,~,~] = fgetMetabolicEnergySmooth2004all(...
                a_col_opt_unsc(count,:)',a_col_opt_unsc(count,:)',...
                full(lMtilde_opt_all),...
                full(vM_opt_all),full(Fce_opt_all)',full(Fpass_opt_all)',...
                MuscleMass.MassM',pctsts,full(Fiso_opt_all)',body_mass,10);
            elseif  strcmp(S.EModel,'Marg1968')
                e_tot_all = fgetMetabolicEnergy_MargariaSmooth(full(Fce_opt_all)',full(vM_opt_all));
            end
        else
            [e_tot_all,~,~,~,~,~] = fgetMetabolicEnergySmooth2004all(...
                a_col_opt_unsc(count,:)',a_col_opt_unsc(count,:)',...
                full(lMtilde_opt_all),...
                full(vM_opt_all),full(Fce_opt_all)',full(Fpass_opt_all)',...
                MuscleMass.MassM',pctsts,full(Fiso_opt_all)',body_mass,10);
        end
        e_tot_opt_all = full(e_tot_all)';
        
        % objective function
        J_opt = J_opt + 1/(dist_trav_opt)*(...
            W.E*B(j+1) * (f_J92exp(e_tot_opt_all,exp_E))/body_mass*h_opt + ...
            W.A*B(j+1) * (f_J92(a_col_opt(count,:)))*h_opt +...
            W.ArmE*B(j+1) * (f_J8(e_a_opt(k,:)))*h_opt +...
            W.Mtp*B(j+1) * (f_J2(e_mtp_opt(k,:)))*h_opt +...
            W.Ak*B(j+1) * (f_J23(qdotdot_col_opt(count,residuals_noarmsi)))*h_opt +...
            W.passMom*B(j+1)* (f_J23(Tau_passj_J(count,:)))*h_opt + ...
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
            (f_J23(Tau_passj_J(count,:)))*h_opt;
        vA_cost = vA_cost + W.u*B(j+1)*...
            (f_J92(vA_opt(k,:)))*h_opt;
        dFTtilde_cost = dFTtilde_cost + W.u*B(j+1)*...
            (f_J92(dFTtilde_col_opt(count,:)))*h_opt;
        QdotdotArm_cost = QdotdotArm_cost + W.u*B(j+1)*...
            (f_J8(qdotdot_col_opt(count,armsi)))*h_opt;
        count = count + 1;
    end
end
J_optf = full(J_opt);           Obj.J = J_optf;
E_costf = full(E_cost);         Obj.E = E_costf;
A_costf = full(A_cost);         Obj.A = A_costf;
Arm_costf = full(Arm_cost);     Obj.Arm = Arm_costf;
Mtp_costf = full(Mtp_cost);     Obj.Mtp = Mtp_costf;
Qdotdot_costf = full(Qdotdot_cost); Obj.qdd = Qdotdot_costf;
Pass_costf = full(Pass_cost);   Obj.Pass = Pass_costf;
vA_costf = full(vA_cost);       Obj.vA = vA_costf;
dFTtilde_costf = full(dFTtilde_cost); Obj.dFTtilde = dFTtilde_costf;
QdotdotArm_costf = full(QdotdotArm_cost); Obj.qdd_arm = QdotdotArm_costf;
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
if ~isfield(S,'Symmetric')
    S.Symmetric = true;
end
if S.Symmetric
    % detect heelstrike
    [IC1i_c,IC1i_s,HS1] = getHeelstrikeSimulation(GRFk_opt,N);
        
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
    dt_exoShift = IC1i_c.*nanmean(diff(tgrid));
    if strcmp(HS1,'l')
        T_exo_GC = T_exo_GC(:,[2 1]);
        dt_exoShift = dt_exoShift - tgrid(end);
    end
    
    % If exoskeleton is implemented as torque actuator
    % evaluate influence of ankle moment and subtalar moment
    if S.ExoBool == 1 && strcmp(ExoImplementation,'TorqueTibiaCalcn')
        % get ID with exoskeleton as percentage of gait cycle
        Ts_opt_Exo = zeros(N*2,size(Qs_opt,2));
        Ts_opt_Exo(1:N-IC1i_c+1,1:nq.all) = Foutk_opt_Exo(IC1i_c:end,1:nq.all);
        Ts_opt_Exo(N-IC1i_c+2:N-IC1i_c+1+N,QsSymA_ptx) = Foutk_opt_Exo(1:end,QsSymB_ptx); %
        Ts_opt_Exo(N-IC1i_c+2:N-IC1i_c+1+N,QsOpp) = -Foutk_opt_Exo(1:end,QsOpp);
        Ts_opt_Exo(N-IC1i_c+2+N:2*N,1:nq.all) = Foutk_opt_Exo(1:IC1i_c-1,1:nq.all);
        % If the first heel strike was on the left foot then we invert so that
        % we always start with the right foot, for analysis purpose
        if strcmp(HS1,'l')
            Ts_opt_Exo(:,QsSymA_ptx) = Ts_opt_Exo(:,QsSymB_ptx);
            Ts_opt_Exo(:,QsOpp) = -Ts_opt_Exo(:,QsOpp);
        end
        Ts_opt_Exo = Ts_opt_Exo./body_mass;
        % compute relative difference
        TExo_Joint = Ts_opt - Ts_opt_Exo;
    end
    
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
    
    
elseif S.Periodic
    % to Do: implement arrange results in gait cycle (should be easy)
    
    % (1) detect heelstrike
    
    % (2) extract states and controls
    
    
end




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
    
    if strcmp(ExoImplementation,'TorqueTibiaCalcn')  || F1.nnz_out == 73
        % compute COP information
        nfr = length(Qs_GC(:,1));
        qdqdd = zeros(nfr,nq.all*2);
        qdqdd(:,1:2:62) = Qs_GC;
        qdqdd(:,2:2:62) = Qdots_GC;
        qdd = Qdotdots_GC;
        qdqdd(:,[1:6 13:end]) = qdqdd(:,[1:6 13:end])*pi./180;
        qdqdd(:,11) = qdqdd(:,11);
        COPR = zeros(nfr,3);    FR = zeros(nfr,3);  MR = zeros(nfr,3);
        COPL = zeros(nfr,3);    FL = zeros(nfr,3);  ML = zeros(nfr,3);
        for ind = 1:nfr
            if F1.nnz_in > nq.all*3
                res = full(F1([qdqdd(ind,:)'; qdd(ind,:)'; ExoZeroT])); % torque L and R exo added
            else
                res = full(F1([qdqdd(ind,:)'; qdd(ind,:)'])); % normal implementation
            end
            % compute the COP position
            FR(ind,:) = res(GRFi.r);
            FL(ind,:) = res(GRFi.l);
            MR(ind,:) = res(GroundT.r);
            ML(ind,:) = res(GroundT.l);
            if abs(FR(ind,2)) > 10
                COPR(ind,:) = [MR(ind,3)./FR(ind,2), 0, -MR(ind,1)./FR(ind,2)];
            end
            if abs(FL(ind,2)) > 10
                COPL(ind,:) = [ML(ind,3)./FL(ind,2), 0, -ML(ind,1)./FL(ind,2)];
            end
        end
        data = [FR COPR FL COPL zeros(nfr,6)];
        dataOut = [data; zeros(size(data))];
        colnames = get_GRFlabels();
        filenameGRF = fullfile(OutFolder,[S.savename '_GRF.mot']);
        time =JointAngleMuscleAct.data(:,1);
        generateMotFile([time dataOut], ['time ' colnames], filenameGRF);
    end
end

% Run bodykinematics
% hard coded for now, this is not a nice implementation
if ~isfield(S,'OsimModel')
    S.OsimModel = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\PredSim_3D\Results\s1_Poggensee_Visual.osim';
end

if exist(S.OsimModel,'file')
   % run body kinematics
   event = [q_opt_GUI_GC(1,1) q_opt_GUI_GC(end,1)];
   [BodyKin] = OpenSim_BodyKinematics(S.OsimModel,pwd,filenameJointAngles,event,true);
else
    BodyKin =[];
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
metab_Etot  = zeros(2*N, NMuscle);
metab_Adot  = zeros(2*N, NMuscle);
metab_Mdot  = zeros(2*N, NMuscle);
metab_Sdot  = zeros(2*N, NMuscle);
metab_Wdot  = zeros(2*N, NMuscle);
FT_opt      = zeros(2*N, NMuscle);
lMT_Vect    = zeros(2*N,92);
vMT_Vect    = zeros(2*N,92);
dM_Vect     = zeros(2*N,92,10);
Fce_opt     = zeros(2*N, NMuscle);
vM_Vect     = zeros(2*N, NMuscle);
Fpass_opt   = zeros(2*N, NMuscle);

metab_Bargh2004  = zeros(2*N,1);  metab_Bargh2004B  = zeros(2*N,1);
metab_Umb2003  = zeros(2*N,1);    metab_Umb2003B  = zeros(2*N,1);
metab_Umb2010  = zeros(2*N,1);    metab_Umb2010B  = zeros(2*N,1);
metab_Uchida2016 = zeros(2*N,1);  metab_Uchida2016B  = zeros(2*N,1);
metab_Umb2010_h1 = zeros(2*N,1);  metab_Umb2010_h1B = zeros(2*N,1);
metab_Umb2010_neg = zeros(2*N,1); metab_Umb2010_negB = zeros(2*N,1);
metab_Marg1968 = zeros(2*N, NMuscle);

for nn = 1:2*N
    % Get muscle-tendon lengths, velocities, moment arms
    % Left leg
    qin_l_opt = Qs_opt_rad(nn,IndexLeft);
    qdotin_l_opt = qdot_opt_GC_rad(nn,IndexLeft);
    [lMTk_l_opt,vMTk_l_opt,dM_l] = f_lMT_vMT_dM(qin_l_opt,qdotin_l_opt);
    % Right leg
    qin_r_opt = Qs_opt_rad(nn,IndexRight);
    qdotin_r_opt = qdot_opt_GC_rad(nn,IndexRight);
    [lMTk_r_opt,vMTk_r_opt,dM_r] = f_lMT_vMT_dM(qin_r_opt,qdotin_r_opt);
    % Both legs
    lMTk_lr_opt     = [lMTk_l_opt([1:43,47:49],1);lMTk_r_opt(1:46,1)];
    vMTk_lr_opt     = [vMTk_l_opt([1:43,47:49],1);vMTk_r_opt(1:46,1)];
    dM_lr_opt       = [dM_l([1:43,47:49],:); dM_r(1:46,:)];
    % force equilibrium
    [~,FT_optt,Fce_optt,Fpass_optt,Fiso_optt] =...
        f_forceEquilibrium_FtildeState_all_tendon(...
        Acts_GC(nn,:)',FTtilde_GC(nn,:)',dFTtilde_GC(nn,:)',full(lMTk_lr_opt),...
        full(vMTk_lr_opt),tensions);
    % fiber kinematics
    [~,lMtilde_opt] = f_FiberLength_TendonForce_tendon(...
        FTtilde_GC(nn,:)',full(lMTk_lr_opt));
    lMtilde_opt_all(nn,:) = full(lMtilde_opt)';
    [vM_opt,vMtilde_opt] = f_FiberVelocity_TendonForce_tendon(FTtilde_GC(nn,:)',...
        dFTtilde_GC(nn,:)',full(lMTk_lr_opt),full(vMTk_lr_opt));
    vMtilde_opt_all(nn,:) = full(vMtilde_opt)';
    
    % Bhargava et al. (2004)
    [energy_total,Adot,Mdot,Sdot,Wdot,eBargh] = ...
        fgetMetabolicEnergySmooth2004all(Acts_GC(nn,:)',...
        Acts_GC(nn,:)',full(lMtilde_opt),full(vM_opt),...
        full(Fce_optt),full(Fpass_optt),MuscleMass.MassM',pctsts,...
        full(Fiso_optt)',body_mass,10);
    
    % Umberger 2003
    vMtildeUmbk_opt = full(vM_opt)./(MTparameters_m(2,:)');
    [eUmb2003,~,~,~,eUmb2003B] = fgetMetabolicEnergySmooth2003all(...
        Acts_GC(nn,:)',Acts_GC(nn,:)',full(lMtilde_opt),...
        vMtildeUmbk_opt,full(vM_opt),full(Fce_optt)',...
        MuscleMass.MassM',pctsts,10,...
        full(Fiso_optt)',body_mass,10);
    
    % Umberger 2010
    [eUmb2010,~,~,~,eUmb2010B] = fgetMetabolicEnergySmooth2010all(...
        Acts_GC(nn,:)',Acts_GC(nn,:)',full(lMtilde_opt),...
        vMtildeUmbk_opt,full(vM_opt),full(Fce_optt)',...
        MuscleMass.MassM',pctsts,10,...
        full(Fiso_optt)',body_mass,10);
    
    % Uchida et al. (2016)
    [eUchida2016,~,~,~,eUchida2016B] = fgetMetabolicEnergySmooth2016all(...
        Acts_GC(nn,:)',Acts_GC(nn,:)',full(lMtilde_opt),...
        vMtildeUmbk_opt,full(vM_opt),full(Fce_optt)',...
        MuscleMass.MassM',pctsts,10,...
        full(Fiso_optt)',body_mass,10);
    
    % Umberger (2010) treating muscle lengthening
    % heat rate as Umberger et al. (2003)
    % vMtilde defined for this model as vM/lMopt
    [eUmb2010_h1,~,~,~,eUmb2010_h1B] = fgetMetabolicEnergySmooth2010all_hl(...
        Acts_GC(nn,:)',Acts_GC(nn,:)',full(lMtilde_opt),...
        vMtildeUmbk_opt,full(vM_opt),full(Fce_optt)',...
        MuscleMass.MassM',pctsts,10,...
        full(Fiso_optt)',body_mass,10);
    
    % Umberger (2010) treating negative mechanical
    % work as Umberger et al. (2003)
    % vMtilde defined for this model as vM/lMopt
    [eUmb2010_neg,~,~,~,eUmb2010_negB] = fgetMetabolicEnergySmooth2010all_neg(...
        Acts_GC(nn,:)',Acts_GC(nn,:)',full(lMtilde_opt),...
        vMtildeUmbk_opt,full(vM_opt),full(Fce_optt)',...
        MuscleMass.MassM',pctsts,10,...
        full(Fiso_optt)',body_mass,10);
    
    % Margaria 1968
    eMarg1968 = getMetabolicEnergy_MargariaSmooth(full(Fce_optt)',full(vM_opt)',1000);
    
    % store results
    e_mo_opt(nn) = full(eBargh)';
    e_mo_optb(nn) = full(eBargh)';
    metab_Etot(nn,:) = full(energy_total)';
    metab_Adot(nn,:) = full(Adot)';
    metab_Mdot(nn,:) = full(Mdot)';
    metab_Sdot(nn,:) = full(Sdot)';
    metab_Wdot(nn,:) = full(Wdot)';
    FT_opt(nn,:)     = full(FT_optt)';
    Fce_opt(nn,:)    = full(Fce_optt)';
    
    lMT_Vect(nn,:)  = full(lMTk_lr_opt);
    vMT_Vect(nn,:)  = full(vMTk_lr_opt);
    dM_Vect(nn,:,:) = full(dM_lr_opt);
    vM_Vect(nn,:)   = full(vM_opt)';
    Fpass_opt(nn,:) = full(Fpass_optt)';
    
    metab_Bargh2004(nn)   = full(sum(energy_total));
    metab_Umb2003(nn)     = full(sum(eUmb2003));
    metab_Umb2010(nn)     = full(sum(eUmb2010));
    metab_Uchida2016(nn)  = full(sum(eUchida2016));
    metab_Umb2010_h1(nn)  = full(sum(eUmb2010_h1));
    metab_Umb2010_neg(nn) = full(sum(eUmb2010_neg));
    metab_Bargh2004B(nn)   = full(eBargh);
    metab_Umb2003B(nn)     = full(eUmb2003B);
    metab_Umb2010B(nn)     = full(eUmb2010B);
    metab_Uchida2016B(nn)  = full(eUchida2016B);
    metab_Umb2010_h1B(nn)  = full(eUmb2010_h1B);
    metab_Umb2010_negB(nn) = full(eUmb2010_negB);
    metab_Marg1968(nn,:)    = full(eMarg1968);
    
end
% Get COT
dist_trav_opt_GC = Qs_opt_rad(end,jointi.pelvis.tx) - ...
    Qs_opt_rad(1,jointi.pelvis.tx); % distance traveled
time_GC = q_opt_GUI_GC(:,1);
e_mo_opt_trb = trapz(time_GC,e_mo_optb);
% Cost of transport: J/kg/m
% Energy model from Bhargava et al. (2004)
COT_GC = e_mo_opt_trb/body_mass/dist_trav_opt_GC;

% COT for all models
COTv.Bargh2004 = trapz(time_GC,metab_Bargh2004)/body_mass/dist_trav_opt_GC;
COTv.Umb2003 = trapz(time_GC,metab_Umb2003)/body_mass/dist_trav_opt_GC;
COTv.Umb2010 = trapz(time_GC,metab_Umb2010)/body_mass/dist_trav_opt_GC;
COTv.Uchida2016 = trapz(time_GC,metab_Uchida2016)/body_mass/dist_trav_opt_GC;
COTv.Umb2010_h1 = trapz(time_GC,metab_Umb2010_h1)/body_mass/dist_trav_opt_GC;
COTv.Umb2010_neg = trapz(time_GC,metab_Umb2010_neg)/body_mass/dist_trav_opt_GC;
COTv.Marg1968 = sum(trapz(time_GC,metab_Marg1968))/body_mass/dist_trav_opt_GC;

% COT for all models (with basal Rate
COTvB.Bargh2004 = trapz(time_GC,metab_Bargh2004B)/body_mass/dist_trav_opt_GC;
COTvB.Umb2003 = trapz(time_GC,metab_Umb2003B)/body_mass/dist_trav_opt_GC;
COTvB.Umb2010 = trapz(time_GC,metab_Umb2010B)/body_mass/dist_trav_opt_GC;
COTvB.Uchida2016 = trapz(time_GC,metab_Uchida2016B)/body_mass/dist_trav_opt_GC;
COTvB.Umb2010_h1 = trapz(time_GC,metab_Umb2010_h1B)/body_mass/dist_trav_opt_GC;
COTvB.Umb2010_neg = trapz(time_GC,metab_Umb2010_negB)/body_mass/dist_trav_opt_GC;
COTvB.Marg1968 = sum(trapz(time_GC,metab_Marg1968))/body_mass/dist_trav_opt_GC;

% Store Energy
EnergyV.Bargh2004       = metab_Bargh2004;
EnergyV.Umb2003         = metab_Umb2003;
EnergyV.Umb2010         = metab_Umb2010;
EnergyV.Uchida2016      = metab_Uchida2016;
EnergyV.Umb2010_h1      = metab_Umb2010_h1;
EnergyV.Umb2010_neg     = metab_Umb2010_neg;
EnergyV.Marg1968        = sum(metab_Marg1968,2);

% Store Energy (with basal rate)
EnergyVB.Bargh2004       = metab_Bargh2004B;
EnergyVB.Umb2003         = metab_Umb2003B;
EnergyVB.Umb2010         = metab_Umb2010B;
EnergyVB.Uchida2016      = metab_Uchida2016B;
EnergyVB.Umb2010_h1      = metab_Umb2010_h1B;
EnergyVB.Umb2010_neg     = metab_Umb2010_negB;
EnergyVB.Marg1968        = sum(metab_Marg1968,2);

%% subtract energy cost of standing in the computations of COT
% Note: this is de default method in pulmonary gas exchange papers.

MetabStandingFile =  fullfile(pathRepo,'StaticStanding','MetabRate_Standing_s1_Poggensee.mat');
if exist(MetabStandingFile,'file')
    Rstanding           = load(MetabStandingFile);
    COTrel.Bargh2004	= trapz(time_GC,(metab_Bargh2004 - Rstanding.Edot.Bargh2004))/body_mass/dist_trav_opt_GC;
    COTrel.Umb2003      = trapz(time_GC,(metab_Umb2003 - Rstanding.Edot.eUmb2003))/body_mass/dist_trav_opt_GC;
    COTrel.Umb2010      = trapz(time_GC,(metab_Umb2010 - Rstanding.Edot.eUmb2010))/body_mass/dist_trav_opt_GC;
    COTrel.Uchida2016   = trapz(time_GC,(metab_Uchida2016 - Rstanding.Edot.eUchida2016))/body_mass/dist_trav_opt_GC;
    COTrel.Umb2010_h1   = trapz(time_GC,(metab_Umb2010_h1 - Rstanding.Edot.eUmb2010_h1))/body_mass/dist_trav_opt_GC;
    COTrel.Umb2010_neg  = trapz(time_GC,(metab_Umb2010_neg - Rstanding.Edot.eUmb2010_neg))/body_mass/dist_trav_opt_GC;
    COTrel.Marg1968     = sum(trapz(time_GC,metab_Marg1968 - Rstanding.Edot.eMarg1968))/body_mass/dist_trav_opt_GC;
else
    COTrel.Bargh2004	= [];
    COTrel.Umb2003      = [];
    COTrel.Umb2010      = [];
    COTrel.Uchida2016   = [];
    COTrel.Umb2010_h1   = [];
    COTrel.Umb2010_neg  = [];
    COTrel.Marg1968     = [];
end

%% Analyse passive exoskeleton support
% Nuckols 2019
if strcmp(ExoImplementation,'Nuckols2019')
    % compute torque of AFO
    GRFrY = GRFs_opt(:,2).*body_weight/100;
    GRFlY = GRFs_opt(:,5).*body_weight/100;
    
    % detect stance phase
    StanceR = tanh(0.1*(GRFrY-30))*0.5+0.5;
    StanceL = tanh(0.1*(GRFlY-30))*0.5+0.5;
    
    % detect angles above threshold
    ql = Qs_GC(:,jointi.ankle.l)*pi/180;
    qr = Qs_GC(:,jointi.ankle.r)*pi/180;
    BoolqL = tanh(200*(ql-S.AFO_q0*pi/180))*0.5+0.5;
    BoolqR = tanh(200*(qr-S.AFO_q0*pi/180))*0.5+0.5;
    TAFO_l = -S.AFO_stiffness.*StanceL.*BoolqL.*(ql-S.AFO_q0.*pi/180);
    TAFO_r = -S.AFO_stiffness.*StanceR.*BoolqR.*(qr-S.AFO_q0.*pi/180);
    
end

%% Mechanical energy analysis
%------------------------------

% Compute the mechanical work done by the model.
% Computations are done here on the interpolated gait cycle data.
% Note that this might introduce some small inconsistensies.

tau = Ts_opt.*body_mass;
qd = Qdots_GC.*pi./180;

% power for each joint
Power = tau .* qd;

% total mechanical power
TotalPower = sum(Power,2);

% mechanical work at joints
t = q_opt_GUI_GC(:,1);
JointWork = trapz(t,Power);
TotalWork = trapz(t,TotalPower);

% muscle power
vM = -vM_Vect;
Fce = Fce_opt;
MusclePower = Fce.*vM;

% Muscle work
MuscleWorkInd = trapz(t,MusclePower);
MuscleWork    = sum(MuscleWorkInd);

% positive and negative muscle power
PosPower = MusclePower;
PosPower(PosPower<0) = 0;
NegPower = MusclePower;
NegPower(NegPower>0) = 0;
MusclePosWork = trapz(t,PosPower);
MuscleNegWork = trapz(t,NegPower);
MusclePosWorkTotal = sum(MusclePosWork);
MuscleNegWorkTotal = sum(MuscleNegWork);

% store outcomes
R.MechE.JointPower = Power;
R.MechE.TotalPower = TotalPower;
R.MechE.JointWork = JointWork;
R.MechE.TotalWork = TotalWork;
R.MechE.MusclePower = MusclePower;
R.MechE.MuscleWork = MuscleWork;
R.MechE.PosPower = PosPower;
R.MechE.NegPower = NegPower;
R.MechE.MusclePosWork = MusclePosWork;
R.MechE.MuscleNegWork = MuscleNegWork;
R.MechE.MusclePosWorkTotal = MusclePosWorkTotal;
R.MechE.MuscleNegWorkTotal = MuscleNegWorkTotal;


%% Save results
% Structure Results_all
R.t_step    = tgrid;
R.tf_step   = tgrid(end);
R.t         = q_opt_GUI_GC(:,1);
R.tend      = q_opt_GUI_GC(end,1) - q_opt_GUI_GC(1,1);
R.Qs        = Qs_GC;
R.Qdots     = Qdots_GC;
R.Qddots    = Qdotdots_GC;
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
% R.ExoControl  = ExoControl;
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
R.dt_exoShift = dt_exoShift;
R.Obj         = Obj;
R.lMT         = lMT_Vect;
R.vMT         = vMT_Vect;
R.dM          = dM_Vect;
R.Muscle.Fce  = Fce_opt;
R.Muscle.vM   = vM_Vect;
R.Muscle.Fpas = Fpass_opt;
R.Muscle.FT   = FT_opt;
R.Muscle.MTparameters = MTparameters_m;
R.COTv        = COTv;
R.Energy      = EnergyV;
R.COTv_basal  = COTvB;
R.Energy_basal= EnergyVB;
R.COTrel      = COTrel;
R.stats       = stats;
R.BodyKin     = BodyKin;

if S.ExoBool == 1 && strcmp(ExoImplementation,'TorqueTibiaCalcn')
    R.TidExo = Ts_opt_Exo.*body_mass;
    R.Exodiff_id = TExo_Joint.*body_mass;
end

if F1.nnz_out == 73
    R.COPL = COPL;
    R.COPR = COPR;
end

% nuckols 2019 results
if strcmp(ExoImplementation,'Nuckols2019')
    R.T_exo       = [TAFO_l TAFO_r];
    R.Nuckols.q0  = S.AFO_q0;
    R.Nuckols.k   = S.AFO_stiffness;
else
    R.Nuckols = [];
end
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
R.colheaders.dM = {'hip flex','hip add','hip rot','knee angle','ankle angle','subtalar angle','mtp angle', 'trunk ext ','trunk bend','trunk rot'};

%% Additional outcomes

% percentage stance and swing phase
[R.Event.Stance, R.Event.Swing, R.Event.DS] = GetPercentageStance(R.GRFs(:,[2 5]).*body_weight/100,30);

% Stepwidth
if isfield(R,'COPL') && isfield(R,'COPR')
    % compute average positin during left stance
    COPR_mean = nanmean(R.COPR(R.GRFs(:,2).*body_weight/100>30,3));
    COPL_mean = nanmean(R.COPL(R.GRFs(:,5).*body_weight/100>30,3));
    % stepwidth
    R.StepWidth_COP = abs(COPR_mean-COPL_mean);
end
%% Save data
% script information
R.info.script = 'f_LoadSim_PoggenSee2020.m';
% Save data
OutFolder = fullfile(pathRepo,'Results',S.ResultsFolder);
FilenameAnalysis = fullfile(OutFolder,[S.savename '_pp.mat']);
save(FilenameAnalysis,'R');

end

