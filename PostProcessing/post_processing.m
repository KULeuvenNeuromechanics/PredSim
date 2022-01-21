function [R]=post_processing(S,R,model_info,f_casadi)
% post_processing.m
%   Function that runs post processing.
%
% INPUT:
%   ResultsFolder
%   * Folder containing results
%   S
%   * Struct containing settings
%   f_casadi
%   * Structu containing casadi functions
%
% OUTPUT:
%   R
%   * Struct containing results

% Original author: Maarten Afschrift
% Original date: ?
%
% Last edit by: Tom Buurke
% Last edit date: 02/12/2021

%% Settings
% Import Casadi
import casadi.*

% Settings
body_mass   = S.mass;
body_weight = S.mass*9.81;
v_tgt       = S.v_tgt;      % average speed
N           = S.N;          % number of mesh intervals
W           = S.W;
exp_E       = S.W.exp_E;    % power metabolic energy
tol_ipopt   = S.tol_ipopt;  % ipopt options

% Model info
nq          = model_info.ExtFunIO.nq;% Number of degrees of freedom

%% Extract necessary parameters from R
% To be written, R structure unclear to me.

%% Load external functions
% The external function performs inverse dynamics through the
% OpenSim/Simbody C++ API. This external function is compiled as a dll from
% which we create a Function instance using CasADi in MATLAB. More details
% about the external function can be found in the documentation.
% Loading external functions.
pathmain = pwd;
[filepath,~,~] =  fileparts(pathmain);
[pathRepo,~,~] = fileparts(filepath);
pathExternalFunctions = [pathRepo,'/ExternalFunctions'];
cd(pathExternalFunctions);
F  = external('F',S.ExternalFunc);
cd(pathmain);

%% load the metalbolic energy equations
PathEnergyEq = fullfile(pathCasADiFunctions,'EnergyModels');
cd(PathEnergyEq);
fgetMetabolicEnergySmooth2003all    = Function.load('fgetMetabolicEnergySmooth2003all');
fgetMetabolicEnergySmooth2010all    = Function.load('fgetMetabolicEnergySmooth2010all');
fgetMetabolicEnergySmooth2016all    = Function.load('fgetMetabolicEnergySmooth2016all');
fgetMetabolicEnergySmooth2010all_hl = Function.load('fgetMetabolicEnergySmooth2010all_hl');
fgetMetabolicEnergySmooth2010all_neg= Function.load('fgetMetabolicEnergySmooth2010all_neg');
cd(pathmain);

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
R.dist_r = sqrt(f_Jnn3(out_res_opt_all(end,calcOrall.r)-...
    out_res_opt_all(1,calcOrall.r)));
% Left leg
R.dist_l = sqrt(f_Jnn3(out_res_opt_all(end,calcOrall.l)-...
    out_res_opt_all(1,calcOrall.l)));
% The total stride length is the sum of the right and left stride
% lengths after a half gait cycle, since we assume symmetry
R.StrideLength_opt = full(dist_r + dist_l);
% The stride width is the medial distance between the calcaneus origins
R.StepWidth_opt = full(abs(out_res_opt_all(:,calcOrall.r(3)) - ...
    out_res_opt_all(:,calcOrall.l(3))));
R.stride_width_mean = mean(StepWidth_opt);

%% Assert average speed
%%% Not sure whether this is used? %%%
dist_trav_opt = q_opt_unsc_all.rad(end,jointi.pelvis.tx) - ...
    q_opt_unsc_all.rad(1,jointi.pelvis.tx); % distance traveled
time_elaps_opt = tf_opt; % time elapsed
vel_aver_opt = dist_trav_opt/time_elaps_opt;
% assert_v_tg should be 0
assert_v_tg = abs(vel_aver_opt-v_tgt);
if assert_v_tg > 1*10^(-tol_ipopt)
    disp('Issue when reconstructing average speed')
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

%% Mechanical energy analysis
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
% Structure R
R.t_step    = tgrid;
R.tf_step   = tgrid(end);
R.t         = q_opt_GUI_GC(:,1);
R.tend      = q_opt_GUI_GC(end,1) - q_opt_GUI_GC(1,1);

R.Ts        = Ts_opt;
R.Tid       = Ts_opt.*body_mass;

R.COT       = COT_GC;
R.StrideLength = StrideLength_opt;
R.StepWidth = stride_width_mean;
R.vMtilde   = vMtilde_opt_all;
R.lMtilde   = lMtilde_opt_all;
R.MetabB.Etot = metab_Etot;
R.MetabB.Adot = metab_Adot;
R.MetabB.Mdot = metab_Mdot;
R.MetabB.Sdot = metab_Sdot;

R.body_mass   = body_mass;
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
% Percentage stance and swing phase
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
OutFolder = fullfile(pathRepo,'Results',S.ResultsFolder);
FilenameAnalysis = fullfile(OutFolder,[S.savename '_pp.mat']);
save(FilenameAnalysis,'R');

end
