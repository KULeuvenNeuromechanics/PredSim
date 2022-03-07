function [] = post_processing(model_info,S,f_casadi)
%%
%% User inputs (typical settings structure)
% settings for optimization
N           = S.solver.N_meshes;          % number of mesh intervals
W           = S.weights;          % weights optimization

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
pathSubjectFolder = S.subject.save_results;
cd(pathSubjectFolder)
F  = external('F',S.ExternalFunc);
cd(pathmain);

%% output folder
OutFolder = S.subject.save_results;

%% Collocation Scheme
% We use a pseudospectral direct collocation method, i.e. we use Lagrange
% polynomials to approximate the state derivatives at the collocation
% points in each mesh interval. We use d=3 collocation points per mesh
% interval and Radau collocation points.
d = 3; % degree of interpolating polynomial
method = 'radau'; % collocation method
[tau_root,C,D,B] = CollocationScheme(d,method);

%% Muscle information
% Muscles from one leg and from the back
muscleNames = model_info.muscle_info.muscle_names;

% Total number of muscles
NMuscle = model_info.muscle_info.NMuscle;
[~,mai] = MomentArmIndices_asym(muscleNames,...
    model_info.muscle_info.polyFit.muscle_spanning_joint_info);
% calculate total number of joints that each muscle crosses (used later)
sumCross = sum(model_info.muscle_info.polyFit.muscle_spanning_joint_info);

% Parameters for activation dynamics
tact = 0.015; % Activation time constant
tdeact = 0.06; % Deactivation time constant

%% Metabolic energy model parameters
% We extract the specific tensions and slow twitch rations.
tensions = model_info.muscle_info.tensions;
pctsts = model_info.muscle_info.pctsts;

%% Function to compute muscle mass
MuscleMass = model_info.muscle_info.muscle_mass;

% Outname = fullfile(OutFolder,[S.subject.name '.mat']);
% load(Outname,'w_opt','stats','setup','Sopt');

%% Stride length and width
% For the stride length we also need the values at the end of the
% interval so N+1 where states but not controls are defined
Xk_Qs_Qdots_opt_all = zeros(N+1,2*size(q_opt_unsc_all.rad,2));
Xk_Qs_Qdots_opt_all(:,1:2:end)  = q_opt_unsc_all.rad;
Xk_Qs_Qdots_opt_all(:,2:2:end)  = qdot_opt_unsc_all.rad;
% We just want to extract the positions of the calcaneus origins so we
% do not really care about Qdotdot that we set to 0
Xk_Qdotdots_opt_all = zeros(N+1,size(q_opt_unsc_all.rad,2));
out_res_opt_all = zeros(N+1,F.nnz_out);
ndof = size(q_opt_unsc_all.rad,2);
for i = 1:N+1
    [res] = F([Xk_Qs_Qdots_opt_all(i,:)';Xk_Qdotdots_opt_all(i,:)']);
    out_res_opt_all(i,:) = full(res);
end
% The stride length is the distance covered by the calcaneus origin
% Right leg
dist_r = sqrt(f_casadi.J_nn_3(out_res_opt_all(end,model_info.ExtFunIO.origin.calcn_r)-...
    out_res_opt_all(1,model_info.ExtFunIO.origin.calcn_r)));
% Left leg
dist_l = sqrt(f_casadi.J_nn_3(out_res_opt_all(end,model_info.ExtFunIO.origin.calcn_l)-...
    out_res_opt_all(1,model_info.ExtFunIO.origin.calcn_l)));
if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    % The total stride length is the sum of the right and left stride
    % lengths after a half gait cycle, since we assume symmetry
    StrideLength_opt = full(dist_r + dist_l);
else
    StrideLength_opt = full(dist_r + dist_l)/2;
end    
% The stride width is the medial distance between the calcaneus origins
StepWidth_opt = full(abs(out_res_opt_all(:,model_info.ExtFunIO.origin.calcn_r(3)) - ...
    out_res_opt_all(:,model_info.ExtFunIO.origin.calcn_l(3))));
stride_width_mean = mean(StepWidth_opt);
stride_width_std = std(StepWidth_opt);

%%
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
Qs_opt_rad(:,model_info.ExtFunIO.jointi.rotations) = Qs_opt_rad(:,model_info.ExtFunIO.jointi.rotations).*pi/180;
qdot_opt_GC_rad = Qdots_GC;
qdot_opt_GC_rad(:,model_info.ExtFunIO.jointi.rotations)= qdot_opt_GC_rad(:,model_info.ExtFunIO.jointi.rotations).*pi/180;
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

for nn = 1:2*N
    % Get muscle-tendon lengths, velocities, moment arms
    % Left leg
    qin_opt = Qs_opt_rad(nn,:);
    qdotin_opt = qdot_opt_GC_rad(nn,:);
    [lMTk_opt,vMTk_opt,dM_opt] = f_casadi.lMT_vMT_dM(qin_opt,qdotin_opt);
    % force equilibrium
    [~,FT_optt,Fce_optt,Fpass_optt,Fiso_optt] =...
        f_casadi.forceEquilibrium_FtildeState_all_tendon(...
        Acts_GC(nn,:)',FTtilde_GC(nn,:)',dFTtilde_GC(nn,:)',full(lMTk_opt),...
        full(vMTk_opt),tensions);
    % fiber kinematics
    [~,lMtilde_opt] = f_casadi.FiberLength_TendonForce_tendon(...
        FTtilde_GC(nn,:)',full(lMTk_opt));
    lMtilde_opt_all(nn,:) = full(lMtilde_opt)';
    [vM_opt,vMtilde_opt] = f_casadi.FiberVelocity_TendonForce_tendon(FTtilde_GC(nn,:)',...
        dFTtilde_GC(nn,:)',full(lMTk_opt),full(vMTk_opt));
    vMtilde_opt_all(nn,:) = full(vMtilde_opt)';
    
    % Bhargava et al. (2004)
    [energy_total,Adot,Mdot,Sdot,Wdot,eBargh] = ...
        f_casadi.getMetabolicEnergySmooth2004all(Acts_GC(nn,:)',...
        Acts_GC(nn,:)',full(lMtilde_opt),full(vM_opt),...
        full(Fce_optt),full(Fpass_optt),MuscleMass',pctsts,...
        full(Fiso_optt)',S.subject.mass,10);
    
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
    
    lMT_Vect(nn,:)  = full(lMTk_opt);
    vMT_Vect(nn,:)  = full(vMTk_opt);
    dM_Vect(nn,:,:) = full(dM_opt);
    vM_Vect(nn,:)   = full(vM_opt)';
    Fpass_opt(nn,:) = full(Fpass_optt)';
    
    metab_Bargh2004(nn)   = full(sum(energy_total));
    metab_Bargh2004B(nn)   = full(eBargh);
end
% Get COT
dist_trav_opt_GC = Qs_opt_rad(end,model_info.ExtFunIO.coordi.pelvis_tx) - ...
    Qs_opt_rad(1,model_info.ExtFunIO.coordi.pelvis_tx); % distance traveled
time_GC = q_opt_GUI_GC(:,1);
e_mo_opt_trb = trapz(time_GC,e_mo_optb);
% Cost of transport: J/kg/m
% Energy model from Bhargava et al. (2004)
COT_GC = e_mo_opt_trb/S.subject.mass/dist_trav_opt_GC;

% COT for all models
COTv.Bargh2004 = trapz(time_GC,metab_Bargh2004)/S.subject.mass/dist_trav_opt_GC;

% COT for all models (with basal Rate
COTvB.Bargh2004 = trapz(time_GC,metab_Bargh2004B)/S.subject.mass/dist_trav_opt_GC;

% Store Energy
EnergyV.Bargh2004       = metab_Bargh2004;

% Store Energy (with basal rate)
EnergyVB.Bargh2004       = metab_Bargh2004B;

%% subtract energy cost of standing in the computations of COT
% Note: this is de default method in pulmonary gas exchange papers.

MetabStandingFile =  fullfile(pathRepo,'StaticStanding','MetabRate_Standing_s1_Poggensee.mat');
if exist(MetabStandingFile,'file')
    Rstanding           = load(MetabStandingFile);
    COTrel.Bargh2004	= trapz(time_GC,(metab_Bargh2004 - Rstanding.Edot.Bargh2004))/S.subject.mass/dist_trav_opt_GC;
else
    COTrel.Bargh2004	= [];
end

%% Mechanical energy analysis
%------------------------------

% Compute the mechanical work done by the model.
% Computations are done here on the interpolated gait cycle data.
% Note that this might introduce some small inconsistensies.

tau = Ts_opt.*S.subject.mass;
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

%% store outcomes
FilenameAnalysis = fullfile(OutFolder,[S.subject.name '_pp.mat']);
laod(FilenameAnalysis,'R');

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
R.a         = Acts_GC;
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
R.body_mass   = S.subject.mass;
R.FT          = FT_opt;
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
R.BodyKin     = BodyKin;

% header information
R.colheaders.joints = joints;
R.colheaders.GRF = {'fore_aft_r','vertical_r',...
    'lateral_r','fore_aft_l','vertical_l','lateral_l'};
R.colheaders.muscles = model_info.muscle_info.muscle_names;
R.colheaders.dM = fields(model_info.ExtFunIO.coordi);

%% Additional outcomes

% percentage stance and swing phase
[R.Event.Stance, R.Event.Swing, R.Event.DS] = GetPercentageStance(R.GRFs(:,[2 5]).*body_weight/100,30);

%% Save data
% Save data
save(FilenameAnalysis,'R');
