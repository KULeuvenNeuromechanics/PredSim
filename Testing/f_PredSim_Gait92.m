function [] = f_PredSim_Gait92(model_info,S,f_casadi)
%%
% TO CHECK
% MTP is not a part of passive torques that go in the cost funciton, Antoine had it in there
% scaling.tauMTP = 100 is used, Antoine had it 30
% Eccentric muscle force was not a part of Antoine's implementation
% Why is muscle torque not a part of this MTP torque path constraint?
% Not implemented minimum foot height during swing, need origin velocity indecies

% %% Adding the casadi path seems to be needed to run processes in batch
% AddCasadiPaths();

% %% Default settings
% S = GetDefaultSettings(S);

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
testing=1;
if testing
F  = external('F',S.ExternalFunc);
else
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

% We use a pseudospectral direct collocation method, i.e. we use Lagrange
% polynomials to approximate the state derivatives at the collocation
% points in each mesh interval. We use d=3 collocation points per mesh
% interval and Radau collocation points.
pathCollocationScheme = [pathRepo,'/CollocationScheme'];
addpath(genpath(pathCollocationScheme));

% Muscle indices for later use
pathmusclemodel = [pathRepo,'/MuscleModel'];
addpath(genpath(pathmusclemodel));

pathMetabolicEnergy = [pathRepo,'/MetabolicEnergy'];
addpath(genpath(pathMetabolicEnergy));

% The initial guess depends on the settings
pathIG = [pathRepo,'/IG'];
addpath(genpath(pathIG));
pathIG = [pathRepo,'/IG'];
addpath(genpath(pathIG));

end

d = 3; % degree of interpolating polynomial
method = 'radau'; % collocation method
[~,C,D,B] = CollocationScheme(d,method);

%% Muscle information
% Muscles from one leg and from the back
muscleNames = model_info.muscle_info.params.names;

% Total number of muscles
NMuscle = size(model_info.muscle_info.params.params,2);
[~,mai] = MomentArmIndices_asym(muscleNames,...
    model_info.muscle_info.polyFit.muscle_spanning_joint_info);
sumCross = sum(model_info.muscle_info.polyFit.muscle_spanning_joint_info);

% Parameters for activation dynamics
tact = 0.015; % Activation time constant
tdeact = 0.06; % Deactivation time constant

%% Metabolic energy model parameters
% We extract the specific tensions and slow twitch rations.
tensions = getSpecificTensions(muscleNames); % CAN ADD TO PREPROCESSING
% tensions = [tension;tension];
pctsts = getSlowTwitchRatios(muscleNames); % CAN ADD TO PREPROCESSING
% pctsts = [pctst;pctst];

% %% CasADi functions
% % We create several CasADi functions for later use
% pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
% PathDefaultFunc = fullfile(pathCasADiFunctions,S.CasadiFunc_Folders);
% f_ArmActivationDynamics = Function.load(fullfile(PathDefaultFunc,'f_ArmActivationDynamics'));
% f_FiberLength_TendonForce_tendon = Function.load(fullfile(PathDefaultFunc,'f_FiberLength_TendonForce_tendon'));
% f_FiberVelocity_TendonForce_tendon = Function.load(fullfile(PathDefaultFunc,'f_FiberVelocity_TendonForce_tendon'));
% f_forceEquilibrium_FtildeState_all_tendon = Function.load(fullfile(PathDefaultFunc,'f_forceEquilibrium_FtildeState_all_tendon'));
% f_J2    = Function.load(fullfile(PathDefaultFunc,'f_J2'));
% f_J23   = Function.load(fullfile(PathDefaultFunc,'f_J23'));
% f_J8    = Function.load(fullfile(PathDefaultFunc,'f_J8'));
% f_J92   = Function.load(fullfile(PathDefaultFunc,'f_J92'));
% f_J92exp = Function.load(fullfile(PathDefaultFunc,'f_J92exp'));
% f_Jnn2  = Function.load(fullfile(PathDefaultFunc,'f_Jnn2'));
% f_lMT_vMT_dM = Function.load(fullfile(PathDefaultFunc,'f_lMT_vMT_dM'));
% f_MtpActivationDynamics = Function.load(fullfile(PathDefaultFunc,'f_MtpActivationDynamics'));
% f_T12 = Function.load(fullfile(PathDefaultFunc,'f_T12'));
% f_T13 = Function.load(fullfile(PathDefaultFunc,'f_T13'));
% f_T27 = Function.load(fullfile(PathDefaultFunc,'f_T27'));
% f_T6 = Function.load(fullfile(PathDefaultFunc,'f_T6'));
% f_AllPassiveTorques = Function.load(fullfile(PathDefaultFunc,'f_AllPassiveTorques'));
% fgetMetabolicEnergySmooth2004all = Function.load(fullfile(PathDefaultFunc,'fgetMetabolicEnergySmooth2004all'));
% 
% PathEnergyEq = fullfile(pathCasADiFunctions,'EnergyModels');
% fgetMetabolicEnergy_MargariaSmooth  = Function.load(fullfile(PathEnergyEq, 'fgetMetabolicEnergy_MargariaSmooth'));

%% Function to compute muscle mass
% BoolRajagopal = 1;
MuscleMass = GetMuscleMass(muscleNames,model_info.muscle_info.params.params);

%% Get bounds and initial guess
if testing
    load('IK_Bounds_Default.mat')
    Qs_walk = getIK_testing(Qsall,model_info);
else
% Kinematics file for bounds -- input arguments
IKfile_bounds = fullfile(pathRepo, S.subject.IG_bounds);

% We extract experimental data to set bounds and initial guesses if needed
joints = fields(model_info.ExtFunIO.coordi)';
% joints = {'pelvis_tilt','pelvis_list','pelvis_rotation','pelvis_tx',...
%     'pelvis_ty','pelvis_tz','hip_flexion_l','hip_adduction_l',...
%     'hip_rotation_l','hip_flexion_r','hip_adduction_r','hip_rotation_r',...
%     'knee_angle_l','knee_angle_r','ankle_angle_l','ankle_angle_r',...
%     'subtalar_angle_l','subtalar_angle_r','mtp_angle_l','mtp_angle_r',...
%     'lumbar_extension','lumbar_bending','lumbar_rotation','arm_flex_l',...
%     'arm_add_l','arm_rot_l','arm_flex_r','arm_add_r','arm_rot_r',...
%     'elbow_flex_l','elbow_flex_r'};
Qs_walk          = getIK(IKfile_bounds,model_info);
end

if S.Bounds_Running 
    [bounds,scaling] = getBounds_all_RunningMaartenB(Qs_walk,model_info,S); % NOTE: used scaling.mtpTau = 100, Antoine used 30
else
    [bounds,scaling] = getBounds_all(Qs_walk,model_info,S); % NOTE: used scaling.mtpTau = 100, Antoine used 30
end

% AdaptBounds.m NO LONGER NEEDED NOW. It was adapting lower bound activation of all
% muscles to a certain value and the bounds on final time. We decided to
% remove lower bound activation of all muscles for now, I implemented the
% bounds on the final time already in the getBounds function
% bounds = AdaptBounds(bounds,S,mai);

% if S.IGsel == 1 % Quasi-random initial guess
%     guess = getGuess_QR_opti_int(N,nq,NMuscle,scaling,S.v_tgt,jointi,d,S.IG_PelvisY);
% elseif S.IGsel == 2 % Data-informed initial guess
%     if S.IGmodeID  < 3 % Data from average walking motion
%         IKfile_guess    = fullfile(pathRepo, S.IKfile_guess);
%         Qs_guess        = getIK(IKfile_guess,joints);
%         time_IC         = [Qs_guess.time(1),Qs_guess.time(end)];
%         guess = getGuess_DI_opti_int_mtp(Qs_guess,nq,N,time_IC,NMuscle,jointi,...
%             scaling,S.v_tgt,d);
%     elseif S.IGmodeID == 3 || S.IGmodeID == 4 % Data from selected motion
%         % Extract joint positions from existing motion (previous results)
%         if S.IGmodeID == 3
%             GuessFolder = fullfile(pathRepo,'Results',S.ResultsF_ig);
%         elseif S.IGmodeID ==4
%             GuessFolder = fullfile(pathRepo,'IG','data');
%         end
%         pathIK      = fullfile(GuessFolder,[S.savename_ig '.mot']);
%         Qs_ig       = getIK(pathIK,joints);
%         % When saving the results, we save a 2 full gait cycle (4*N) so here we
%         % only select 1:N to have half a gait cycle
%         nfr = length(Qs_ig.allfilt(:,1));
%         frSel = round(nfr./4);
%         Qs_ig_sel.allfilt   = Qs_ig.allfilt(1:frSel,:);
%         Qs_ig_sel.time      = Qs_ig.time(1:frSel,:);
%         Qs_ig_sel.colheaders = Qs_ig.colheaders;
%         time_IC = [Qs_ig_sel.time(1),Qs_ig_sel.time(end)];
%         guess = getGuess_DI_opti_int_mtp(Qs_ig_sel,nq,N,time_IC,NMuscle,jointi,scaling,S.v_tgt,d);
%     end
% end

if strcmp(S.subject.IG_selection,'quasi-random')
    guess = getGuess_QR_opti(N,scaling,model_info,S,d);
else
    IKfile_guess    = fullfile(pathRepo, S.subject.IG_selection);
    Qs_guess        = getIK(IKfile_guess,model_info);
    time_IC         = [Qs_guess.time(1),Qs_guess.time(end)];
    guess = getGuess_DI_opti(Qs_guess,N,time_IC,scaling,S,d,model_info);
end

% adapt guess so that it fits within the bounds
guess = AdaptGuess_UserInput(guess,bounds,S);

nq = model_info.ExtFunIO.nq;
visualizebounds

bounds.tf.lower = guess.tf;
bounds.tf.upper = guess.tf;
% %% exoskeleton torques
% % function to get exoskeleton torques at mesh points
% if ~S.ExoDesignBool
%     ExoVect = GetExoTorques(S,pathRepo,N);
% end


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
lboundsQsk(model_info.ExtFunIO.coordi.pelvis_tx,1) = ...
    bounds.QsQdots_0.lower(2*model_info.ExtFunIO.coordi.pelvis_tx-1);
uboundsQsk = bounds.QsQdots.upper(1:2:end)'*ones(1,N+1);
uboundsQsk(model_info.ExtFunIO.coordi.pelvis_tx,1) = ...
    bounds.QsQdots_0.upper(2*model_info.ExtFunIO.coordi.pelvis_tx-1);
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
if strcmp(S.subject.mtp_type,'active')
a_mtp = opti.variable(nq.mtp,N+1);
opti.subject_to(bounds.a_mtp.lower'*ones(1,N+1) < a_mtp < ...
    bounds.a_mtp.upper'*ones(1,N+1));
opti.set_initial(a_mtp, guess.a_mtp');
% Mtp activations at collocation points
a_mtp_col = opti.variable(nq.mtp,d*N);
end
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Design optimal exoskeleton support
% if S.OptTexo_Ankle.Bool
%     ExoVect = opti.variable(2, N);
%     opti.subject_to(S.OptTexo_Ankle.Tbound(1)< ExoVect < S.OptTexo_Ankle.Tbound(2));
%     opti.set_initial(ExoVect,zeros(2,N));
% elseif S.OptTexo_AnkleKneeHip.Bool
%     ExoVect = opti.variable(6, N);
%     opti.subject_to(S.OptTexo_AnkleKneeHip.Tbound_Ankle(1) < ...
%         ExoVect(1:2,:) < S.OptTexo_AnkleKneeHip.Tbound_Ankle(2));
%     opti.subject_to(S.OptTexo_AnkleKneeHip.Tbound_Knee(1) < ...
%         ExoVect(3:4,:) < S.OptTexo_AnkleKneeHip.Tbound_Knee(2));
%     opti.subject_to(S.OptTexo_AnkleKneeHip.Tbound_Hip(1) < ...
%         ExoVect(5:6,:) < S.OptTexo_AnkleKneeHip.Tbound_Hip(2));
% end

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

% % define the exoskeleton assistive torque
% nExoDofs = length(ExoVect(:,1));
% Texok   = MX.sym('Texo',nExoDofs,1); % joint moments for the exoskeleton

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
%     qinj    = Qskj_nsc(model_info.ExtFunIO.jointi.legs_torso, j+1);
%     qdotinj = Qdotskj_nsc(model_info.ExtFunIO.jointi.legs_torso, j+1);
    qinj    = Qskj_nsc(:, j+1);
    qdotinj = Qdotskj_nsc(:, j+1);
    [lMTj,vMTj,MAj] =  f_casadi.lMT_vMT_dM(qinj',qdotinj');
    for i=1:model_info.ExtFunIO.nq.legs_torso
        MA_j.(model_info.ExtFunIO.jointi.names.legs_torso{i}) = ...
            MAj(mai(model_info.ExtFunIO.jointi.legs_torso(i)).mus',i);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get muscle-tendon forces and derive Hill-equilibrium
    [Hilldiffj,FTj,Fcej,Fpassj,Fisoj] = ...
        f_casadi.forceEquilibrium_FtildeState_all_tendon(akj(:,j+1),...
        FTtildekj_nsc(:,j+1),dFTtildej_nsc(:,j),...
        lMTj,vMTj,tensions);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get metabolic energy rate if in the cost function
    % Get muscle fiber lengths
    [~,lMtildej] = f_casadi.FiberLength_TendonForce_tendon(...
        FTtildekj_nsc(:,j+1),lMTj);
    % Get muscle fiber velocities
    [vMj,~] = f_casadi.FiberVelocity_TendonForce_tendon(...
        FTtildekj_nsc(:,j+1),dFTtildej_nsc(:,j),...
        lMTj,vMTj);
    if strcmp(S.EModel,'Bhargava2004')
        % Get metabolic energy rate Bhargava et al. (2004)
        [e_totj,~,~,~,~,~] = f_casadi.getMetabolicEnergySmooth2004all(...
            akj(:,j+1),akj(:,j+1),lMtildej,vMj,Fcej,Fpassj,...
            MuscleMass',pctsts,Fisoj,S.subject.mass,10);
    elseif  strcmp(S.EModel,'Marg1968')
        e_totj = f_casadi.getMetabolicEnergy_MargariaSmooth(Fcej,vMj);
    else
        error('No energy model selected');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get passive joint torques
    Tau_passj_all = f_casadi.AllPassiveTorques(Qskj_nsc(:,j+1),Qdotskj_nsc(:,j+1));
    Tau_passj_arm = Tau_passj_all(:,(model_info.ExtFunIO.nq.muscleActuated+1):(model_info.ExtFunIO.nq.muscleActuated+model_info.ExtFunIO.nq.arms));
    Tau_passj_mtp = Tau_passj_all(:,(model_info.ExtFunIO.nq.muscleActuated+model_info.ExtFunIO.nq.arms+1):end);
    Tau_passj_noMTP = Tau_passj_all(:,1:(model_info.ExtFunIO.nq.muscleActuated+model_info.ExtFunIO.nq.arms));% MTP is not a part of this, Antoine had it in there
    Tau_passj_muscleActuated = Tau_passj_all(:,1:model_info.ExtFunIO.nq.muscleActuated);
    
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
    da_adtj = f_casadi.ArmActivationDynamics(e_ak,a_akj(:,j+1)');
    eq_constr{end+1} = (h*da_adtj - a_ap)./scaling.a_a;
    % Mtp activation dynamics (explicit formulation)
    da_mtpdtj = f_casadi.MtpActivationDynamics(e_mtpk,a_mtpkj(:,j+1)');
    eq_constr{end+1} = (h*da_mtpdtj - a_mtpp);
    % Add contribution to the quadrature function
%     J = J + 1*(...
%         W.E*B(j+1)      *(f_J92exp(e_totj,W.exp_E))/S.mass*h + ...
%         W.A*B(j+1)      *(f_J92(akj(:,j+1)'))*h + ...
%         W.ArmE*B(j+1)   *(f_J8(e_ak))*h +...
%         W.Mtp*B(j+1)    *(f_J2(e_mtpk))*h +...
%         W.Ak*B(j+1)     *(f_J23(Aj(coord_noarmsi,j)))*h + ...
%         W.passMom*B(j+1)*(f_J23(Tau_passj_J))*h + ...
%         W.u*B(j+1)      *(f_J92(vAk))*h + ...
%         W.u*B(j+1)      *(f_J92(dFTtildej(:,j)))*h + ...
%         W.u*B(j+1)      *(f_J8(Aj(armsi,j)))*h);

    J = J + 1*(...
        W.E*B(j+1)          *(f_casadi.J_N_muscles_exp(e_totj,W.E_exp))/S.subject.mass*h + ...
        W.a*B(j+1)          *(f_casadi.J_N_muscles(akj(:,j+1)'))*h + ...
        W.e_arm*B(j+1)      *(f_casadi.J_arms_dof(e_ak))*h +...
        W.e_mtp*B(j+1)      *(f_casadi.J_2(e_mtpk))*h +...
        W.q_dotdot*B(j+1)   *(f_casadi.J_noarms_dof(Aj(model_info.ExtFunIO.jointi.noarmsi,j)))*h + ...
        W.pass_torq*B(j+1)  *(f_casadi.J_pass_noMTP_dof(Tau_passj_noMTP))*h + ... % MTP is not a part of this, Antoine had it in there
        W.slack_ctrl*B(j+1) *(f_casadi.J_N_muscles(vAk))*h + ...
        W.slack_ctrl*B(j+1) *(f_casadi.J_N_muscles(dFTtildej(:,j)))*h + ...
        W.slack_ctrl*B(j+1) *(f_casadi.J_arms_dof(Aj(model_info.ExtFunIO.jointi.armsi,j)))*h);

    if W.EccF ~= 0 % was not a part of Antoine's implementation
        % select eccentric muscle force
        b = 100;
        ivMPos = (tanh(b*(vMj-0.02))*0.5+0.5);      % negative power = 0, positive power = positive power
        Fecc = FTj.*ivMPos;
        % add Fecc squared to the objective function
        J = J + 1*(...
            W.EccF*B(j+1)    *(f_casadi.J_N_muscles_exp(Fecc,2))*h );
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Call external function (run inverse dynamics)
    if F.nnz_in == model_info.ExtFunIO.nq.all*3
        % no exoskeleton torque as input in passive simulations
        [Tj] = F([QsQdotskj_nsc(:,j+1);Aj_nsc(:,j)]);    % left and right leg exoskeleton torques as inputs as well.
    elseif F.nnz_in > model_info.ExtFunIO.nq.all*3
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
    eq_constr{end+1} = Tj(model_info.ExtFunIO.jointi.ground_pelvis,1);
    % Muscle-driven joint torques for the lower limbs and the trunk
    for i=1:model_info.ExtFunIO.nq.muscleActuated
        Ft.(model_info.ExtFunIO.jointi.names.muscleActuated{i}) = FTj(mai(model_info.ExtFunIO.jointi.muscleActuated(i)).mus',1);
        T.(model_info.ExtFunIO.jointi.names.muscleActuated{i}) = ...
            f_casadi.(['musc_cross_' num2str(sumCross(model_info.ExtFunIO.jointi.muscleActuated(i)))])...
            (MA_j.(model_info.ExtFunIO.jointi.names.muscleActuated{i}),...
            Ft.(model_info.ExtFunIO.jointi.names.muscleActuated{i}));
        eq_constr{end+1} = Tj(model_info.ExtFunIO.jointi.muscleActuated(i),1) - ...
            (T.(model_info.ExtFunIO.jointi.names.muscleActuated{i}) + ...
            Tau_passj_muscleActuated(:,i));
    end
    % Arms
    eq_constr{end+1} = Tj(model_info.ExtFunIO.jointi.armsi,1)/scaling.ArmTau - (a_akj(:,j+1) + ...
        (Tau_passj_arm')/scaling.ArmTau);
    % Mtp
    eq_constr{end+1} = Tj(model_info.ExtFunIO.jointi.mtpi,1)/scaling.MtpTau - (a_mtpkj(:,j+1) + ...
        (Tau_passj_mtp')/scaling.MtpTau); % Why is muscle torque not a part of this?
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
    Qconstr = f_casadi.J_nn_2(Tj(model_info.ExtFunIO.origin.calcn_r([1 3]),1) - ...
        Tj(model_info.ExtFunIO.origin.calcn_l([1 3]),1));
    ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) = Qconstr;
    ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 3;
    ctIneq = ctIneq + length(Qconstr);
    % Constraint to prevent the arms to penetrate the skeleton
    % Origins femurs and ipsilateral hands (transv plane) at minimum
    % 18 cm from each other.
    Qconstr = f_casadi.J_nn_2(Tj(model_info.ExtFunIO.origin.femur_r([1 3]),1) - ...
        Tj(model_info.ExtFunIO.origin.hand_r([1 3]),1));
    ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) = Qconstr;
    ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 4;
    ctIneq = ctIneq + length(Qconstr);
    Qconstr = f_casadi.J_nn_2(Tj(model_info.ExtFunIO.origin.femur_l([1 3]),1) - ...
        Tj(model_info.ExtFunIO.origin.hand_l([1 3]),1));
    ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) = Qconstr;
    ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 4;
    ctIneq = ctIneq + length(Qconstr);
    % Origins tibia (transv plane) at minimum 11 cm from each other.
    Qconstr = f_casadi.J_nn_2(Tj(model_info.ExtFunIO.origin.tibia_r([1 3]),1) - ...
        Tj(model_info.ExtFunIO.origin.tibia_l([1 3]),1));
    ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) = Qconstr;
    ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 5;
    ctIneq = ctIneq + length(Qconstr);
    % Origins toes (transv plane) at minimum 10 cm from each other.
    Qconstr = f_casadi.J_nn_2(Tj(model_info.ExtFunIO.origin.toes_r([1 3]),1) - ...
        Tj(model_info.ExtFunIO.origin.toes_l([1 3]),1));
    ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) = Qconstr;
    ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 6;
    ctIneq = ctIneq + length(Qconstr);
    
%     % Minimum foot height during swing
% CAN BE UNCOMMENTED, NEED INDECIES OF ORIGIN VELOCITIES
%     if isfield(S.Constr,'ImposeFootHeight') && S.Constr.ImposeFootHeight
%         Vx = Tj([toesOr_dot.r(1) toesOr_dot.l(1) calcOr_dot.r(1) calcOr_dot.l(1)]);
%         Ry = Tj([toesOry.r toesOry.l calcOry.r calcOry.l]);
%         vx_tr = S.Constr.FootHeight.Vx_Treshold;
%         tr_ry = S.Constr.FootHeight.Ry_Treshold;
%         b = S.Constr.FootHeight.b;
%         vx_bool = 0.5.*tanh(b*(Vx-vx_tr))+0.5;
%         ry_ad = -Ry+tr_ry;
%         Qconstr =  ry_ad.*vx_bool;
%         ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) =Qconstr;
%         ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 101;
%         ctIneq = ctIneq + length(Qconstr);
%     end
%     % inequality constraint on exoskeleton moments
%     if S.OptTexo_Ankle.Bool
%         Qconstr = Texok(1).*Qskj_nsc(IO.jointi.ankle.l,j+1);
%         ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) = Qconstr;
%         ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 7;
%         ctIneq = ctIneq + length(Qconstr);
%         Qconstr = Texok(2).*Qskj_nsc(IO.jointi.ankle.r,j+1);
%         ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) = Qconstr;
%         ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 7;
%         ctIneq = ctIneq + length(Qconstr);
%     elseif S.OptTexo_AnkleKneeHip.Bool
%         iVect = [IO.jointi.ankle.l IO.jointi.ankle.r IO.jointi.knee.l IO.jointi.knee.r , ...
%             IO.jointi.hip_flex.l IO.jointi.hip_flex.r];
%         for iv = 1:6
%             Qconstr = Texok(iv).*Qskj_nsc(iVect(iv),j+1);
%             ineq_constr(ctIneq:ctIneq+length(Qconstr)-1) = Qconstr;
%             ineq_constr_index(ctIneq:ctIneq+length(Qconstr)-1) = 7+iv-1;
%             ctIneq = ctIneq + length(Qconstr);
%         end        
%     end
end % End loop over collocation points
%%
eq_constrV = vertcat(eq_constr{:});

% Casadi function to get constraints and objective
% f_coll = Function('f_coll',{tfk,ak,aj,FTtildek,FTtildej,Qsk,Qsj,Qdotsk,...
%     Qdotsj,a_ak,a_aj,a_mtpk,a_mtpj,vAk,e_ak,e_mtpk,dFTtildej,Aj,Texok},...
%     {eq_constrV,ineq_constr,J});
f_coll = Function('f_coll',{tfk,ak,aj,FTtildek,FTtildej,Qsk,Qsj,Qdotsk,...
    Qdotsj,a_ak,a_aj,a_mtpk,a_mtpj,vAk,e_ak,e_mtpk,dFTtildej,Aj},...
    {eq_constrV,ineq_constr,J});
% assign NLP problem to multiple cores
f_coll_map = f_coll.map(N,S.solver.parallel_mode,S.solver.N_threads);
if strcmp(S.subject.mtp_type,'active')
% [coll_eq_constr, coll_ineq_constr, Jall] = f_coll_map(tf,...
%     a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, Qs(:,1:end-1), ...
%     Qs_col, Qdots(:,1:end-1), Qdots_col, a_a(:,1:end-1), a_a_col, ...
%     a_mtp(:,1:end-1), a_mtp_col, vA, e_a, e_mtp, dFTtilde_col, A_col,ExoVect);
[coll_eq_constr, coll_ineq_constr, Jall] = f_coll_map(tf,...
    a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, Qs(:,1:end-1), ...
    Qs_col, Qdots(:,1:end-1), Qdots_col, a_a(:,1:end-1), a_a_col, ...
    a_mtp(:,1:end-1), a_mtp_col, vA, e_a, e_mtp, dFTtilde_col, A_col);
else
[coll_eq_constr, coll_ineq_constr, Jall] = f_coll_map(tf,...
    a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, Qs(:,1:end-1), ...
    Qs_col, Qdots(:,1:end-1), Qdots_col, a_a(:,1:end-1), a_a_col, ...
    vA, e_a, e_mtp, dFTtilde_col, A_col);
end
% equality constrains
opti.subject_to(coll_eq_constr == 0);

% inequality constraints (logical indexing not possible in MX arrays)
cSel = coll_ineq_constr(find(ineq_constr_index == 1),:); % activation dyanmics
opti.subject_to(cSel(:)  >= 0);
cSel = coll_ineq_constr(find(ineq_constr_index == 2),:); % deactivation dyanmics
opti.subject_to(cSel(:)  <= 1/tact);
cSel = coll_ineq_constr(find(ineq_constr_index == 3),:); % origin calcaneus
opti.subject_to(S.bounds.calcn_dist.lower.^2 < cSel(:) < 4);
cSel = coll_ineq_constr(find(ineq_constr_index == 4),:); % arms
opti.subject_to(0.0324 < cSel(:) < 4);
cSel = coll_ineq_constr(find(ineq_constr_index == 5),:); % origin tibia minimum x cm away from each other
opti.subject_to(S.bounds.tibia_dist.lower.^2 < cSel(:) < 4);
cSel = coll_ineq_constr(find(ineq_constr_index == 6),:); % origins toes minimum x cm away from each other
opti.subject_to(S.bounds.toes_dist.lower.^2 < cSel(:) < 4);
% if S.OptTexo_Ankle.Bool         % bound on exoskeleton power
%     cSel = coll_ineq_constr(find(ineq_constr_index == 7),:); % exoskeleton power
%     opti.subject_to(S.OptTexo_Ankle.Pbound(1) < cSel(:) < S.OptTexo_Ankle.Pbound(2));
% elseif S.OptTexo_AnkleKneeHip.Bool
%     lb_ExoPower = [S.OptTexo_AnkleKneeHip.Pbound_Ankle(1)*ones(1,2) ,...
%         S.OptTexo_AnkleKneeHip.Pbound_Knee(1)*ones(1,2), ...
%         S.OptTexo_AnkleKneeHip.Pbound_Hip(1)*ones(1,2)];
%     ub_ExoPower = [S.OptTexo_AnkleKneeHip.Pbound_Ankle(2)*ones(1,2) ,...
%         S.OptTexo_AnkleKneeHip.Pbound_Knee(2)*ones(1,2), ...
%         S.OptTexo_AnkleKneeHip.Pbound_Hip(2)*ones(1,2)];
%     for iv = 1:6
%         cSel = coll_ineq_constr(find(ineq_constr_index == 7+iv-1),:); % exoskeleton power
%         opti.subject_to(lb_ExoPower(iv) < cSel(:) < ub_ExoPower(iv));
%     end
% end

% % inequality constraint minimum foot height during swing
% CAN BE UNCOMMENTED, NEED INDECIES OF ORIGIN VELOCITIES
% if isfield(S.Constr,'ImposeFootHeight') && S.Constr.ImposeFootHeight
%     cSel = coll_ineq_constr(find(ineq_constr_index == 101),:); % origins toes minimum x cm away from each other
%     opti.subject_to(cSel(:) < 0.01);
% end

% Loop over mesh points
for k=1:N
    % Variables within current mesh interval
    % States
    akj = [a(:,k), a_col(:,(k-1)*d+1:k*d)];
    FTtildekj = [FTtilde(:,k), FTtilde_col(:,(k-1)*d+1:k*d)];
    Qskj = [Qs(:,k), Qs_col(:,(k-1)*d+1:k*d)];
    Qdotskj = [Qdots(:,k), Qdots_col(:,(k-1)*d+1:k*d)];
    a_akj = [a_a(:,k), a_a_col(:,(k-1)*d+1:k*d)];
    if strcmp(S.subject.mtp_type,'active')
        a_mtpkj = [a_mtp(:,k), a_mtp_col(:,(k-1)*d+1:k*d)];
        opti.subject_to(a_mtp(:,k+1) == a_mtpkj*D);
    end
    % Add equality constraints (next interval starts with end values of
    % states from previous interval)
    opti.subject_to(a(:,k+1) == akj*D);
    opti.subject_to(FTtilde(:,k+1) == FTtildekj*D); % scaled
    opti.subject_to(Qs(:,k+1) == Qskj*D); % scaled
    opti.subject_to(Qdots(:,k+1) == Qdotskj*D); % scaled
    opti.subject_to(a_a(:,k+1) == a_akj*D);
end % End loop over mesh points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional path constraints
if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    % Periodicity of the states (or rather LR symmetry -half gait cycle)
    % Qs and Qdots
    opti.subject_to(Qs(model_info.ExtFunIO.symQs.QsInvA,end) - Qs(model_info.ExtFunIO.symQs.QsInvB,1) == 0);
    opti.subject_to(Qdots(model_info.ExtFunIO.symQs.QdotsInvA,end) - Qdots(model_info.ExtFunIO.symQs.QdotsInvB,1) == 0);
    opti.subject_to(Qs(model_info.ExtFunIO.symQs.orderQsOpp,end) + Qs(model_info.ExtFunIO.symQs.orderQsOpp,1) == 0);
    opti.subject_to(Qdots(model_info.ExtFunIO.symQs.orderQsOpp,end) + Qdots(model_info.ExtFunIO.symQs.orderQsOpp,1) == 0);
    % Muscle activations
    orderMusInv = [NMuscle/2+1:NMuscle,1:NMuscle/2];
    opti.subject_to(a(:,end) - a(orderMusInv,1) == 0);
    % Muscle-tendon forces
    opti.subject_to(FTtilde(:,end) - FTtilde(orderMusInv,1) == 0);
%     % Arm activations
%     orderArmInv = [IO.jointi.sh_flex.r:IO.jointi.sh_rot.r,...
%         IO.jointi.sh_flex.l:IO.jointi.sh_rot.l,...
%         IO.jointi.elb.r:IO.jointi.elb.r,...
%         IO.jointi.elb.l:IO.jointi.elb.l]-IO.jointi.sh_flex.l+1;
    opti.subject_to(a_a(:,end) - a_a(model_info.ExtFunIO.symQs.orderArmInv,1) == 0);
    % Mtp activations
%     orderMtpInv = [model_info.ExtFunIO.coordi.mtp_angle_r,model_info.ExtFunIO.coordi.mtp_angle_l];
%     orderMtpInv = orderMtpInv-min(orderMtpInv)+1;
    orderMtpInv = [2 1];
    if strcmp(S.subject.mtp_type,'active')
        opti.subject_to(a_mtp(:,end) - a_mtp(orderMtpInv,1) == 0);
    end
else
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
    if strcmp(S.subject.mtp_type,'active')
        opti.subject_to(a_mtp(:,end) - a_mtp(:,1) == 0);
    end
end
% Average speed
% Provide expression for the distance traveled
Qs_nsc = Qs.*(scaling.QsQdots(1:2:end)'*ones(1,N+1));
dist_trav_tot = Qs_nsc(model_info.ExtFunIO.coordi.pelvis_tx,end) - ...
    Qs_nsc(model_info.ExtFunIO.coordi.pelvis_tx,1);
vel_aver_tot = dist_trav_tot/tf;
opti.subject_to(vel_aver_tot - S.subject.v_pelvis_x_trgt == 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale cost function
Jall_sc = sum(Jall)/dist_trav_tot;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create NLP solver
opti.minimize(Jall_sc);
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.mu_strategy           = 'adaptive';
options.ipopt.max_iter              = S.solver.max_iter;
options.ipopt.linear_solver         = S.solver.linear_solver;
options.ipopt.tol                   = 1*10^(-S.solver.tol_ipopt);
options.ipopt.constr_viol_tol       = 1*10^(-S.solver.tol_ipopt);
opti.solver('ipopt', options);
% Create and save diary
if testing
    OutFolder = S.subject.save_folder;
else
OutFolder = fullfile(pathRepo,'Results',S.ResultsFolder);
end
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

