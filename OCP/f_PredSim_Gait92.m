function [] = f_PredSim_Gait92(model_info,S,f_casadi)
%%
% TO CHECK
% MTP is not a part of passive torques that go in the cost funciton, Antoine had it in there
% Not implemented minimum foot height during swing, need origin velocity indecies
% Notes:
% Need External function, IG file in Subjects/<subject_name> folder
% Need collocationscheme.m and other files/functions of MuscleModel and
% MetabolicEnergy in VariousFunctions folder

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
pathSubjectFolder = [pathRepo,'/Subjects/' S.subject.name];
cd(pathSubjectFolder)
F  = external('F',S.ExternalFunc);
cd(pathmain);

% We use a pseudospectral direct collocation method, i.e. we use Lagrange
% polynomials to approximate the state derivatives at the collocation
% points in each mesh interval. We use d=3 collocation points per mesh
% interval and Radau collocation points.
d = 3; % degree of interpolating polynomial
method = 'radau'; % collocation method
[~,C,D,B] = CollocationScheme(d,method);

%% Muscle information
% Muscles from one leg and from the back
muscleNames = model_info.muscle_info.muscle_names;

% Total number of muscles
NMuscle = size(model_info.muscle_info.params,2);
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

%% Get bounds and initial guess
% Kinematics file for bounds -- input arguments
IKfile_bounds = fullfile(pathSubjectFolder, S.subject.IG_bounds);

% We extract experimental data to set bounds and initial guesses if needed
Qs_walk          = getIK(IKfile_bounds,model_info);

if S.Bounds_Running
    % NOTE: Running simulations needs to be tested
    [bounds,scaling] = getBounds_all_Running(Qs_walk,model_info,S);
else
    [bounds,scaling] = getBounds_all(Qs_walk,model_info,S);
end

if strcmp(S.subject.IG_selection,'quasi-random')
    guess = getGuess_QR_opti(N,scaling,model_info,S,d);
else
    IKfile_guess    = fullfile(pathRepo, S.subject.IG_selection);
    Qs_guess        = getIK(IKfile_guess,model_info);
    time_IC         = [Qs_guess.time(1),Qs_guess.time(end)];
    guess = getGuess_DI_opti(Qs_guess,N,time_IC,scaling,S,d,model_info);
end

% adapt guess so that it fits within the bounds
[guess,bounds] = AdaptGuess_UserInput(guess,bounds,S);

nq = model_info.ExtFunIO.jointi.nq;

if (S.visualize_IG_bounds)
    visualizebounds
end

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
% Mtp activations at mesh points
if strcmp(S.subject.mtp_type,'active')
    a_mtp = opti.variable(nq.mtp,N+1);
    opti.subject_to(bounds.a_mtp.lower'*ones(1,N+1) < a_mtp < ...
        bounds.a_mtp.upper'*ones(1,N+1));
    opti.set_initial(a_mtp, guess.a_mtp');
    % Mtp activations at collocation points
    a_mtp_col = opti.variable(nq.mtp,d*N);
    opti.subject_to(bounds.a_mtp.lower'*ones(1,d*N) < a_mtp_col < ...
        bounds.a_mtp.upper'*ones(1,d*N));
    opti.set_initial(a_mtp_col, guess.a_mtp_col');
end
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
if strcmp(S.subject.mtp_type,'active')
    e_mtp = opti.variable(nq.mtp, N);
    opti.subject_to(bounds.e_mtp.lower'*ones(1,N) < e_mtp < ...
        bounds.e_mtp.upper'*ones(1,N));
    opti.set_initial(e_mtp, guess.e_mtp');
end
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

%% OCP: collocation equations
% Define CasADi variables for static parameters
tfk         = MX.sym('tfk'); % MX variable for final time
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
% Define CasADi variables for controls
vAk     = MX.sym('vAk',NMuscle);
e_ak    = MX.sym('e_ak',nq.arms);

if strcmp(S.subject.mtp_type,'active')
    a_mtpk      = MX.sym('a_mtpk',nq.mtp);
    a_mtpj      = MX.sym('a_mtpkmesh',nq.mtp,d);
    a_mtpkj     = [a_mtpk a_mtpj];
    % Define CasADi variables for MTP controls
    e_mtpk  = MX.sym('e_mtpk',nq.mtp);
end

% Define CasADi variables for "slack" controls
dFTtildej   = MX.sym('dFTtildej',NMuscle,d);
Aj          = MX.sym('Aj',nq.all,d);
J           = 0; % Initialize cost function
eq_constr   = {}; % Initialize equality constraint vector
% ineq_constr = MX(567,1); % Initialise inequality constraint vector
% ineq_constr_index = nan(567); % keeps track of the index of the inequality constraint
ineq_constr = MX(1000,1); % Initialise inequality constraint vector
ineq_constr_index = nan(1000,1); % keeps track of the index of the inequality constraint
ctIneq      = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time step
h = tfk/N;
% Loop over collocation points
for j=1:d
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unscale variables
    Qskj_nsc = Qskj.*(scaling.QsQdots(1:2:end)'*ones(1,size(Qskj,2)));
    Qdotskj_nsc = Qdotskj.*(scaling.QsQdots(2:2:end)'*ones(1,size(Qdotskj,2)));
    FTtildekj_nsc = FTtildekj.*(scaling.FTtilde'*ones(1,size(FTtildekj,2)));
    dFTtildej_nsc = dFTtildej.*scaling.dFTtilde;
    Aj_nsc = Aj.*(scaling.Qdotdots'*ones(1,size(Aj,2)));
    vAk_nsc = vAk.*scaling.vA;
    
    QsQdotskj_nsc = MX(nq.all*2, d+1);
    QsQdotskj_nsc(1:2:end,:) = Qskj_nsc;
    QsQdotskj_nsc(2:2:end,:) = Qdotskj_nsc;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get muscle-tendon lengths, velocities, and moment arms
    qinj    = Qskj_nsc(:, j+1);
    qdotinj = Qdotskj_nsc(:, j+1);
    [lMTj,vMTj,MAj] =  f_casadi.lMT_vMT_dM(qinj',qdotinj');
    % Derive the moment arms of all the muscles crossing each joint
    for i=1:model_info.ExtFunIO.jointi.nq.legs_torso
        MA_j.(model_info.ExtFunIO.jointi.names.legs_torso{i}) = ...
            MAj(mai(model_info.ExtFunIO.jointi.legs_torso(i)).mus',...
            model_info.ExtFunIO.jointi.legs_torso(i));
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
    % In the f_casadi.AllPassiveTorques function, the order in which the
    % passive torques are exported is muscleActuated joints, then arms, then mtp
    Tau_passj_arm = Tau_passj_all((model_info.ExtFunIO.jointi.nq.muscleActuated+1):(model_info.ExtFunIO.jointi.nq.muscleActuated+model_info.ExtFunIO.jointi.nq.arms));
    Tau_passj_mtp = Tau_passj_all((model_info.ExtFunIO.jointi.nq.muscleActuated+model_info.ExtFunIO.jointi.nq.arms+1):end);
    Tau_passj_noMTP = Tau_passj_all(1:(model_info.ExtFunIO.jointi.nq.muscleActuated+model_info.ExtFunIO.jointi.nq.arms));% MTP is not a part of this, Antoine had it in there
    Tau_passj_muscleActuated = Tau_passj_all(1:model_info.ExtFunIO.jointi.nq.muscleActuated);
    
    % Expression for the state derivatives at the collocation points
    Qsp_nsc      = Qskj_nsc*C(:,j+1);
    Qdotsp_nsc   = Qdotskj_nsc*C(:,j+1);
    FTtildep_nsc = FTtildekj_nsc*C(:,j+1);
    ap           = akj*C(:,j+1);
    a_ap         = a_akj*C(:,j+1);
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
    
    % Add contribution to the cost function
    J = J + 1*(...
        W.E*B(j+1)          *(f_casadi.J_N_muscles_exp(e_totj,W.E_exp))/S.subject.mass*h + ...
        W.a*B(j+1)          *(f_casadi.J_N_muscles(akj(:,j+1)'))*h + ...
        W.e_arm*B(j+1)      *(f_casadi.J_arms_dof(e_ak))*h +...
        W.q_dotdot*B(j+1)   *(f_casadi.J_noarms_dof(Aj(model_info.ExtFunIO.jointi.noarmsi,j)))*h + ...
        W.pass_torq*B(j+1)  *(f_casadi.J_pass_noMTP_dof(Tau_passj_noMTP))*h + ... % MTP is not a part of this, Antoine had it in there
        W.slack_ctrl*B(j+1) *(f_casadi.J_N_muscles(vAk))*h + ...
        W.slack_ctrl*B(j+1) *(f_casadi.J_N_muscles(dFTtildej(:,j)))*h + ...
        W.slack_ctrl*B(j+1) *(f_casadi.J_arms_dof(Aj(model_info.ExtFunIO.jointi.armsi,j)))*h);
    
    if strcmp(S.subject.mtp_type,'active')
        a_mtpp       = a_mtpkj*C(:,j+1);
        % Mtp activation dynamics (explicit formulation)
        da_mtpdtj = f_casadi.MtpActivationDynamics(e_mtpk,a_mtpkj(:,j+1)');
        eq_constr{end+1} = (h*da_mtpdtj - a_mtpp);
        J = J + 1*(...
            W.e_mtp*B(j+1)      *(f_casadi.J_2(e_mtpk))*h);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Call external function (run inverse dynamics)
    [Tj] = F([QsQdotskj_nsc(:,j+1);Aj_nsc(:,j)]);
    
    % note that this has to be -Texo, since a positive torque around
    % the z-axis of the tibia results in a plantarflexion moment. (and
    % plantarflexion is defined as negative in opensim and in the
    % exo torque vector.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add path constraints
    % Null pelvis residuals
    eq_constr{end+1} = Tj(model_info.ExtFunIO.jointi.ground_pelvis,1);
    % Muscle-driven joint torques for the lower limbs and the trunk
    % Casadi functions are created and called based on the number of
    % muscles each joint crosses
    for i=1:model_info.ExtFunIO.jointi.nq.muscleActuated
        Ft.(model_info.ExtFunIO.jointi.names.muscleActuated{i}) = FTj(mai(model_info.ExtFunIO.jointi.muscleActuated(i)).mus',1);
        T.(model_info.ExtFunIO.jointi.names.muscleActuated{i}) = ...
            f_casadi.(['musc_cross_' num2str(sumCross(model_info.ExtFunIO.jointi.muscleActuated(i)))])...
            (MA_j.(model_info.ExtFunIO.jointi.names.muscleActuated{i}),...
            Ft.(model_info.ExtFunIO.jointi.names.muscleActuated{i}));
        eq_constr{end+1} = Tj(model_info.ExtFunIO.jointi.muscleActuated(i),1) - ...
            (T.(model_info.ExtFunIO.jointi.names.muscleActuated{i}) + ...
            Tau_passj_muscleActuated(i));
    end
    % Arms
    eq_constr{end+1} = Tj(model_info.ExtFunIO.jointi.armsi,1)/scaling.ArmTau - (a_akj(:,j+1) + ...
        (Tau_passj_arm)/scaling.ArmTau);
    % Mtp
    if strcmp(S.subject.mtp_type,'active')
        eq_constr{end+1} = Tj(model_info.ExtFunIO.jointi.mtpi,1)/scaling.MtpTau - (a_mtpkj(:,j+1) + ...
            (Tau_passj_mtp)/scaling.MtpTau);
    else
        eq_constr{end+1} = Tj(model_info.ExtFunIO.jointi.mtpi,1)/scaling.MtpTau - ...
            ((Tau_passj_mtp)/scaling.MtpTau);
    end        
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
end % End loop over collocation points

eq_constrV = vertcat(eq_constr{:});

if strcmp(S.subject.mtp_type,'active')
    % Casadi function to get constraints and objective
    f_coll = Function('f_coll',{tfk,ak,aj,FTtildek,FTtildej,Qsk,Qsj,Qdotsk,...
        Qdotsj,a_ak,a_aj,a_mtpk,a_mtpj,vAk,e_ak,e_mtpk,dFTtildej,Aj},...
        {eq_constrV,ineq_constr,J});
    % assign NLP problem to multiple cores
    f_coll_map = f_coll.map(N,S.solver.parallel_mode,S.solver.N_threads);
    [coll_eq_constr, coll_ineq_constr, Jall] = f_coll_map(tf,...
        a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, Qs(:,1:end-1), ...
        Qs_col, Qdots(:,1:end-1), Qdots_col, a_a(:,1:end-1), a_a_col, ...
        a_mtp(:,1:end-1), a_mtp_col, vA, e_a, e_mtp, dFTtilde_col, A_col);
else
    % Casadi function to get constraints and objective
    f_coll = Function('f_coll',{tfk,ak,aj,FTtildek,FTtildej,Qsk,Qsj,Qdotsk,...
        Qdotsj,a_ak,a_aj,vAk,e_ak,dFTtildej,Aj},...
        {eq_constrV,ineq_constr,J});
    % assign NLP problem to multiple cores
    f_coll_map = f_coll.map(N,S.solver.parallel_mode,S.solver.N_threads);
    [coll_eq_constr, coll_ineq_constr, Jall] = f_coll_map(tf,...
        a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, Qs(:,1:end-1), ...
        Qs_col, Qdots(:,1:end-1), Qdots_col, a_a(:,1:end-1), a_a_col, ...
        vA, e_a, dFTtilde_col, A_col);
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
cSel = coll_ineq_constr(find(ineq_constr_index == 4),:); % arms from hips
opti.subject_to(0.0324 < cSel(:) < 4);
cSel = coll_ineq_constr(find(ineq_constr_index == 5),:); % origin tibia minimum x cm away from each other
opti.subject_to(S.bounds.tibia_dist.lower.^2 < cSel(:) < 4);
cSel = coll_ineq_constr(find(ineq_constr_index == 6),:); % origins toes minimum x cm away from each other
opti.subject_to(S.bounds.toes_dist.lower.^2 < cSel(:) < 4);

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
    opti.subject_to(a(model_info.ExtFunIO.symQs.orderMus,end) - a(model_info.ExtFunIO.symQs.orderMusInv,1) == 0);
    % Muscle-tendon forces
    opti.subject_to(FTtilde(model_info.ExtFunIO.symQs.orderMus,end) - FTtilde(model_info.ExtFunIO.symQs.orderMusInv,1) == 0);
    % Arm activations
    opti.subject_to(a_a(model_info.ExtFunIO.symQs.orderArm,end) - a_a(model_info.ExtFunIO.symQs.orderArmInv,1) == 0);
    % Mtp activations
    if strcmp(S.subject.mtp_type,'active')
        orderMtpInv = [2 1]; % There are only 2 mtps
        opti.subject_to(a_mtp(:,end) - a_mtp(orderMtpInv,1) == 0);
    end
else
    opti.subject_to(Qs(:,end) - Qs(:,1) == 0);
    opti.subject_to(Qdots(:,end) - Qdots(:,1) == 0);
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
%%
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
setup.tolerance.ipopt = S.solver.tol_ipopt;
setup.bounds = bounds;
setup.scaling = scaling;
setup.guess = guess;

%% Save the results

Outname = fullfile(OutFolder,[S.savename '.mat']);
Sopt = S;
save(Outname,'w_opt','stats','setup','Sopt');
end

