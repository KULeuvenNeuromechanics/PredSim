function [R] = f_LoadSim_Rajagopal_Tracking(ResultsFolder,loadname)


%% path information
if strcmp(loadname(end-3:end),'.mat')
    loadname = loadname(1:end-4);
end

pathmain        = mfilename('fullpath');
[filepath,~,~]  =  fileparts(pathmain);
[pathRepo,~,~]  = fileparts(filepath);
OutFolder       = fullfile(pathRepo,'Results',ResultsFolder);
Outname         = fullfile(OutFolder,[loadname '.mat']);
load(Outname,'w_opt','stats','Sopt','setup');
S = Sopt;

body_mass = S.mass;
S.Model.Rajagopal = true;
%% User inputs (typical settings structure)
% load default CasadiFunctions

% flow control
writeIKmotion   = 1; % set to 1 to write .mot file

% settings for optimization
v_tgt       = S.v_tgt;      % average speed
N           = S.N;          % number of mesh intervals
W           = S.W;
exp_E       = S.W.exp_E;    % power metabolic energy
coCont      = S.coCont;     % co-contraction identifier

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
nq.leg      = 7; % #joints needed for polynomials

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
% if isfield(S,'ExternalFunc2')   
%    F1 = external('F',S.ExternalFunc2);
% elseif F.nnz_in ==  nq.all*3
%    F1 = external('F', 'Browning_2008_pp.dll');
% elseif F.nnz_in ==  nq.all*3 + 2   
%    F1 = external('F','SimExo_3D_ExportAll.dll');
%    disp('Selected SimExo_3D_ExportAll.dll as external function because S.ExternalFunc2 was not specified');
% end
F1 = external('F','TrackSim_Subject1_pp.dll');
cd(pathmain);


%% Test the type of simulation
ExoImplementation = 'IdealAnkle'; % (1) ankle actuation (2) passive afo, (3) torque actuator

if strcmp(ExtF,'PredSim_3D_GRF.dll') && isfield(S,'AFO_stiffness')
    % AFO modelled as in Nuckols 2019: Ultrasound imaging links soleus
    % muscle neuromechanics and energetics during human walking with
    % elastic ankle exoskeletons
    disp('Post processing of passive AFO as in Nuckols - ultrasound paper');
    ExoImplementation = 'Nuckols2019';
elseif F.nnz_in == nq.all*3+2
    % After the first simulations report (May 12 2020), we changed the implemenation of the exoskeleton assistance
    % to a torque actuator at the calcaneus and tibia. We also changed from two .dll files (one for optimization and one
    % for post processing) to one .dll file for everything. 
    ExoImplementation = 'TorqueTibiaCalcn';
else
    disp('Default processing with active exoskeleton (Poggensee paper)')
end




%% Second, origins bodies.
% Calcaneus
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
muscleNames = {'addbrev_r','addlong_r','addmagDist_r','addmagIsch_r','addmagMid_r','addmagProx_r',...
    'bflh_r','bfsh_r','edl_r','ehl_r','fdl_r','fhl_r','gaslat_r','gasmed_r','glmax1_r','glmax2_r',...
    'glmax3_r','glmed1_r','glmed2_r','glmed3_r','glmin1_r','glmin2_r','glmin3_r','grac_r','iliacus_r',...
    'perbrev_r','perlong_r','piri_r','psoas_r','recfem_r','sart_r','semimem_r','semiten_r','soleus_r',...
    'tfl_r','tibant_r','tibpost_r','vasint_r','vaslat_r','vasmed_r'};

% Muscle indices for later use
pathmusclemodel = [pathRepo,'/MuscleModel/',subject];
musi = 1:40;
NMuscle = length(muscleNames)*2;
ExtPoly = '';
load([pathmusclemodel,'/MTparameters_',subject, ExtPoly, '.mat']);
MTparameters_m = [MTparameters(:,musi),MTparameters(:,musi)];

% Muscle-tendon parameters. Row 1: maximal isometric forces; Row 2: optimal
% fiber lengths; Row 3: tendon slack lengths; Row 4: optimal pennation
% angles; Row 5: maximal contraction velocities
pathpolynomial = fullfile(pathRepo,'Polynomials',S.subject);
addpath(genpath(pathpolynomial));
tl = load([pathpolynomial,'/muscle_spanning_joint_INFO_',subject,ExtPoly, '.mat']);
[~,mai] = MomentArmIndices(muscleNames(1:end),...
    tl.muscle_spanning_joint_INFO(1:end,:));

% Parameters for activation dynamics
tact = 0.015; % Activation time constant
tdeact = 0.06; % Deactivation time constant

%% Metabolic energy model parameters
% We extract the specific tensions and slow twitch rations.
pathMetabolicEnergy = [pathRepo,'/MetabolicEnergy'];
addpath(genpath(pathMetabolicEnergy));
% (1:end-3), since we do not want to count twice the back muscles
tension = getSpecificTensions(muscleNames);
tensions = [tension;tension];
% (1:end-3), since we do not want to count twice the back muscles
pctst = getSlowTwitchRatios(muscleNames);
pctsts = [pctst;pctst];

%% CasADi functions
pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
PathDefaultFunc = fullfile(pathCasADiFunctions,S.CasadiFunc_Folders);
f_FiberLength_TendonForce_tendon = Function.load(fullfile(PathDefaultFunc,'f_FiberLength_TendonForce_tendon'));
f_FiberVelocity_TendonForce_tendon = Function.load(fullfile(PathDefaultFunc,'f_FiberVelocity_TendonForce_tendon'));
f_forceEquilibrium_FtildeState_all_tendon = Function.load(fullfile(PathDefaultFunc,'f_forceEquilibrium_FtildeState_all_tendon'));
f_J2    = Function.load(fullfile(PathDefaultFunc,'f_J2'));
f_J3    = Function.load(fullfile(PathDefaultFunc,'f_J3'));
f_J23   = Function.load(fullfile(PathDefaultFunc,'f_J23'));
f_J8    = Function.load(fullfile(PathDefaultFunc,'f_J8'));
f_J80   = Function.load(fullfile(PathDefaultFunc,'f_J80'));
f_J80exp = Function.load(fullfile(PathDefaultFunc,'f_J80exp'));
f_Jnn3  = Function.load(fullfile(PathDefaultFunc,'f_Jnn3'));
f_lMT_vMT_dM = Function.load(fullfile(PathDefaultFunc,'f_lMT_vMT_dM'));
f_AllPassiveTorques = Function.load(fullfile(PathDefaultFunc,'f_AllPassiveTorques'));
fgetMetabolicEnergySmooth2004all = Function.load(fullfile(PathDefaultFunc,'fgetMetabolicEnergySmooth2004all'));

%% file with mass of muscles
MassFile = fullfile(PathDefaultFunc,'MassM.mat');
if exist(MassFile,'file')
   MuscleMass = load(MassFile);
else
    MassFile = fullfile(pathCasADiFunctions,'MassM.mat');
    MuscleMass =load(MassFile);
end
%% load the metalbolic energy equations
PathDefaultFunc = fullfile(pathCasADiFunctions,'EnergyModels_Rajagopal');
cd(PathDefaultFunc);
fgetMetabolicEnergySmooth2003all    = Function.load('fgetMetabolicEnergySmooth2003all');
fgetMetabolicEnergySmooth2010all    = Function.load('fgetMetabolicEnergySmooth2010all');
fgetMetabolicEnergySmooth2016all    = Function.load('fgetMetabolicEnergySmooth2016all');
fgetMetabolicEnergySmooth2010all_hl = Function.load('fgetMetabolicEnergySmooth2010all_hl');
fgetMetabolicEnergySmooth2010all_neg= Function.load('fgetMetabolicEnergySmooth2010all_neg');
fgetMetabolicEnergy_MargariaSmooth  = Function.load('fgetMetabolicEnergy_MargariaSmooth');
cd(pathmain);


%% Experimental data
% We extract experimental data to set bounds and initial guesses if needed
joints = {'pelvis_tilt','pelvis_list','pelvis_rotation','pelvis_tx',...
    'pelvis_ty','pelvis_tz','hip_flexion_l','hip_adduction_l',...
    'hip_rotation_l','hip_flexion_r','hip_adduction_r','hip_rotation_r',...
    'knee_angle_l','knee_angle_r','ankle_angle_l','ankle_angle_r',...
    'subtalar_angle_l','subtalar_angle_r','mtp_angle_l','mtp_angle_r',...
    'lumbar_extension','lumbar_bending','lumbar_rotation','arm_flex_l',...
    'arm_add_l','arm_rot_l','arm_flex_r','arm_add_r','arm_rot_r',...
    'elbow_flex_l','elbow_flex_r'};

%% Bounds
scaling = setup.scaling;

%% Index helpers
% get help indexes for left and right leg and for symmetry constraint
[IndexLeft,IndexRight] = GetIndexHelper(S,jointi);

%% Unpack vector with optimization results

% Note: time is not an optimization variable in the tracking simulations
starti = 1;
BoolTrack  = 1;
if ~isfield(setup,'Bools') || setup.Bools.Tracking == false
    NParameters = 1;
    tf_opt = w_opt(1:NParameters);
    starti = NParameters+1;
    BoolTrack = 0;
else
   tf_opt = setup.dt;
   %tf_opt = 0.1;
end
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

a_lumbar_opt = reshape(w_opt(starti:starti+nq.trunk*(N+1)-1),nq.trunk,N+1)';
starti = starti + nq.trunk*(N+1);
a_lumbar_col_opt = reshape(w_opt(starti:starti+nq.trunk*(d*N)-1),nq.trunk,d*N)';
starti = starti + nq.trunk*(d*N);

vA_opt = reshape(w_opt(starti:starti+NMuscle*N-1),NMuscle,N)';
starti = starti + NMuscle*N;
e_a_opt = reshape(w_opt(starti:starti+nq.arms*N-1),nq.arms,N)';
starti = starti + nq.arms*N;
e_mtp_opt = reshape(w_opt(starti:starti+nq.mtp*N-1),nq.mtp,N)';
starti = starti + nq.mtp*N;

e_lumbar_opt = reshape(w_opt(starti:starti+nq.trunk*N-1),nq.trunk,N)';
starti = starti + nq.trunk*N;

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

a_lumbar_mesh_col_opt=zeros(N*(d+1)+1,nq.trunk);
a_lumbar_mesh_col_opt(1:(d+1):end,:)= a_lumbar_opt;

for k=1:N
    rangei = k*(d+1)-(d-1):k*(d+1);
    rangebi = (k-1)*d+1:k*d;
    a_mesh_col_opt(rangei,:) = a_col_opt(rangebi,:);
    FTtilde_mesh_col_opt(rangei,:) = FTtilde_col_opt(rangebi,:);
    Qs_mesh_col_opt(rangei,:) = Qs_col_opt(rangebi,:);
    Qdots_mesh_col_opt(rangei,:) = Qdots_col_opt(rangebi,:);
    a_a_mesh_col_opt(rangei,:) = a_a_col_opt(rangebi,:);
    a_mtp_mesh_col_opt(rangei,:) = a_mtp_col_opt(rangebi,:);
    a_lumbar_mesh_col_opt(rangei,:) = a_lumbar_col_opt(rangebi,:);
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

for i = 1:N
    % ID moments    
    if F1.nnz_in == nq.all*3 + 2
        if S.ExoBool == 1 && strcmp(ExoImplementation,'TorqueTibiaCalcn')
            [res2] = F1([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)'; -ExoVect(:,i)]);
            [res2_or] = F([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)'; -ExoVect(:,i)]); % ext func used in optimization
        end
        [res] = F1([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)'; 0; 0]);
        [res_or] = F([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)'; 0; 0]); % ext func used in optimization
    else
        [res] = F1([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)']);
        [res_or] = F([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)']); % ext func used in optimization
    end
    Foutk_opt(i,:) = full(res);
    Foutk_opt(i,1:nq.all) = full(res_or(1:nq.all)); % extract ID moments from external function used in the optimization
    if S.ExoBool == 1 && strcmp(ExoImplementation,'TorqueTibiaCalcn')
        Foutk_opt_Exo(i,:) = full(res2);
        Foutk_opt_Exo(i,1:nq.all) = full(res2_or(1:nq.all));
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
    if F1.nnz_in == nq.all*3 + 2
        if S.ExoBool == 1 && strcmp(ExoImplementation,'TorqueTibiaCalcn')
            [res2] = F1([Xj_Qs_Qdots_opt(i,:)';Xj_Qdotdots_opt(i,:)'; -ExoVect(:,iMesh)]); 
            [res2_or] = F([Xj_Qs_Qdots_opt(i,:)';Xj_Qdotdots_opt(i,:)'; -ExoVect(:,iMesh)]); 
        end
        [res] = F1([Xj_Qs_Qdots_opt(i,:)';Xj_Qdotdots_opt(i,:)'; 0; 0]);
        [res_or] = F([Xj_Qs_Qdots_opt(i,:)';Xj_Qdotdots_opt(i,:)'; 0; 0]);
    else
        [res] = F1([Xj_Qs_Qdots_opt(i,:)';Xj_Qdotdots_opt(i,:)']);
        [res_or] = F([Xj_Qs_Qdots_opt(i,:)';Xj_Qdotdots_opt(i,:)']);
    end    
    Foutj_opt(i,:) = full(res);
    Foutj_opt(i,1:nq.all) = full(res_or(1:nq.all)); % extract ID moments from external function used in the optimization
    if S.ExoBool == 1 && strcmp(ExoImplementation,'TorqueTibiaCalcn')
        Foutj_opt_Exo(i,:) = full(res2);
        Foutj_opt_Exo(i,1:nq.all) = full(res2(1:nq.all));
    end
    % passive torques
    Tau_passj_opt_all(i,:) = full(f_AllPassiveTorques(q_col_opt_unsc.rad(i,:),qdot_col_opt_unsc.rad(i,:)));
end

%% Create .mot file for OpenSim GUI
if writeIKmotion
    pathOpenSim = [pathRepo,'/OpenSim'];
    addpath(genpath(pathOpenSim));
    labels = {'time','pelvis_tilt','pelvis_list','pelvis_rotation','pelvis_tx',...
    'pelvis_ty','pelvis_tz','hip_flexion_l','hip_adduction_l',...
    'hip_rotation_l','hip_flexion_r','hip_adduction_r','hip_rotation_r',...
    'knee_angle_l','knee_angle_r','ankle_angle_l','ankle_angle_r',...
    'subtalar_angle_l','subtalar_angle_r','mtp_angle_l','mtp_angle_r',...
    'lumbar_extension','lumbar_bending','lumbar_rotation','arm_flex_l',...
    'arm_add_l','arm_rot_l','arm_flex_r','arm_add_r','arm_rot_r',...
    'elbow_flex_l','elbow_flex_r'};
    % Joint angles
    t = tgrid';
    q.data = [tgrid' q_opt_unsc_all.deg];
    q.labels = labels;
    OutFolder = fullfile(pathRepo,'Results',S.ResultsFolder);
    filenameJointAngles = fullfile(OutFolder,[S.savename '.mot']);
    write_motionFile(q, filenameJointAngles);
end


%% Metabolic cost of transport for a gait cycle
Qs_opt_rad = q_opt_unsc_all.rad;
qdot_opt_GC_rad = qdot_opt_unsc_all

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
lMT_Vect    = zeros(2*N,80);
vMT_Vect    = zeros(2*N,80);
dM_Vect     = zeros(2*N,80,7);
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
    lMTk_lr_opt     = [lMTk_l_opt(1:40,1);lMTk_r_opt(1:40,1)];
    vMTk_lr_opt     = [vMTk_l_opt(1:40,1);vMTk_r_opt(1:40,1)];
    dM_lr_opt       = [dM_l(1:40,:); dM_r(1:40,:)];
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
    eMarg1968 = fgetMetabolicEnergy_MargariaSmooth(full(Fce_optt)',full(vM_opt)');    
    
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
EnergyV.Marg1968        = metab_Marg1968;

% Store Energy (with basal rate)
EnergyVB.Bargh2004       = metab_Bargh2004B;
EnergyVB.Umb2003         = metab_Umb2003B;
EnergyVB.Umb2010         = metab_Umb2010B;
EnergyVB.Uchida2016      = metab_Uchida2016B;
EnergyVB.Umb2010_h1      = metab_Umb2010_h1B;
EnergyVB.Umb2010_neg     = metab_Umb2010_negB;
EnergyVB.Marg1968        = metab_Marg1968;

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
R.dt_exoShift = dt_exoShift;
R.Obj         = Obj;
R.lMT         = lMT_Vect;
R.vMT         = vMT_Vect;
R.dM          = dM_Vect;
R.Muscle.Fce  = Fce_opt;
R.Muscle.vM   = vM_Vect;
R.Muscle.Fpas = Fpass_opt;
R.Muscle.FT   = FT_opt;
R.COTv        = COTv;
R.Energy      = EnergyV;
R.COTv_basal  = COTvB;
R.Energy_basal= EnergyVB;
R.COTrel      = COTrel;

if strcmp(ExoImplementation,'TorqueTibiaCalcn')
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

