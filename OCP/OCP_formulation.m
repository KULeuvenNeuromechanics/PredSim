function [] = OCP_formulation(model_info,S,f_casadi)
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
pathSubjectFolder = S.subject.save_results;
cd(pathSubjectFolder)
F  = external('F',S.ExternalFunc);
cd(pathmain);

%% Output folder
OutFolder = S.subject.save_results;
if ~isfolder(OutFolder)
    mkdir(OutFolder);
end

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

%% Get bounds and initial guess
% Kinematics file for bounds -- input arguments
IKfile_bounds = fullfile(pathSubjectFolder, S.subject.IG_bounds);

% We extract experimental data to set bounds and initial guesses if needed
Qs_walk          = getIK(IKfile_bounds,model_info);

if S.misc.running
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

if (S.misc.visualize_IG_bounds)
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
lboundsQsk = bounds.Qs.lower'*ones(1,N+1);
lboundsQsk(model_info.ExtFunIO.coordi.pelvis_tx,1) = ...
    bounds.Qs_0.lower(model_info.ExtFunIO.coordi.pelvis_tx);
uboundsQsk = bounds.Qs.upper'*ones(1,N+1);
uboundsQsk(model_info.ExtFunIO.coordi.pelvis_tx,1) = ...
    bounds.Qs_0.upper(model_info.ExtFunIO.coordi.pelvis_tx);
opti.subject_to(lboundsQsk < Qs < uboundsQsk);
opti.set_initial(Qs, guess.Qs');
% Qs at collocation points
Qs_col = opti.variable(nq.all,d*N);
opti.subject_to(bounds.Qs.lower'*ones(1,d*N) < Qs_col < ...
    bounds.Qs.upper'*ones(1,d*N));
opti.set_initial(Qs_col, guess.Qs_col');
% Qdots at mesh points
Qdots = opti.variable(nq.all,N+1);
opti.subject_to(bounds.Qdots.lower'*ones(1,N+1) < Qdots < ...
    bounds.Qdots.upper'*ones(1,N+1));
opti.set_initial(Qdots, guess.Qdots');
% Qdots at collocation points
Qdots_col = opti.variable(nq.all,d*N);
opti.subject_to(bounds.Qdots.lower'*ones(1,d*N) < Qdots_col < ...
    bounds.Qdots.upper'*ones(1,d*N));
opti.set_initial(Qdots_col, guess.Qdots_col');
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
    Qskj_nsc = Qskj.*(scaling.Qs'*ones(1,size(Qskj,2)));
    Qdotskj_nsc = Qdotskj.*(scaling.Qdots'*ones(1,size(Qdotskj,2)));
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
    for i=1:nq.legs_torso
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
    if strcmp(S.metabolicE.model,'Bhargava2004')
        % Get metabolic energy rate Bhargava et al. (2004)
        [e_totj,~,~,~,~,~] = f_casadi.getMetabolicEnergySmooth2004all(...
            akj(:,j+1),akj(:,j+1),lMtildej,vMj,Fcej,Fpassj,...
            MuscleMass',pctsts,Fisoj,S.subject.mass,10);
    else
        error('No energy model selected');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get passive joint torques
    Tau_passj_all = f_casadi.AllPassiveTorques(Qskj_nsc(:,j+1),Qdotskj_nsc(:,j+1));
    % In the f_casadi.AllPassiveTorques function, the order in which the
    % passive torques are exported is muscleActuated joints, then arms, then mtp
    Tau_passj_arm = Tau_passj_all((nq.muscleActuated+1):(nq.muscleActuated+nq.arms));
    Tau_passj_mtp = Tau_passj_all((nq.muscleActuated+nq.arms+1):end);
    Tau_passj_noMTP = Tau_passj_all(1:(nq.muscleActuated+nq.arms));% MTP is not a part of this, Antoine had it in there
    Tau_passj_muscleActuated = Tau_passj_all(1:nq.muscleActuated);
    
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
    eq_constr{end+1} = (h*qdotj_nsc - Qsp_nsc)./scaling.Qs';
    eq_constr{end+1} = (h*Aj_nsc(:,j) - Qdotsp_nsc)./scaling.Qdots';
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
            W.e_mtp*B(j+1)  *(f_casadi.J_2(e_mtpk))*h);
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
    for i=1:nq.muscleActuated
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
Qs_nsc = Qs.*(scaling.Qs'*ones(1,N+1));
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
Outname = fullfile(OutFolder,[S.subject.name '_log.txt']);
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
Outname = fullfile(OutFolder,[S.subject.name '.mat']);
Sopt = S;
save(Outname,'w_opt','stats','setup','Sopt');


%% Post processing
%% Essential post processing
body_weight = S.subject.mass*9.81;
QsSymA = model_info.ExtFunIO.symQs.QsInvA;
QsSymB = model_info.ExtFunIO.symQs.QsInvB;
QsOpp = model_info.ExtFunIO.symQs.orderQsOpp;
QsSymA_ptx = model_info.ExtFunIO.symQs.QdotsInvA;
QsSymB_ptx = model_info.ExtFunIO.symQs.QdotsInvB;

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
if strcmp(S.subject.mtp_type,'active')
    a_mtp_opt = reshape(w_opt(starti:starti+nq.mtp*(N+1)-1),nq.mtp,N+1)';
    starti = starti + nq.mtp*(N+1);
    a_mtp_col_opt = reshape(w_opt(starti:starti+nq.mtp*(d*N)-1),nq.mtp,d*N)';
    starti = starti + nq.mtp*(d*N);
    a_mtp_mesh_col_opt=zeros(N*(d+1)+1,nq.mtp);
    a_mtp_mesh_col_opt(1:(d+1):end,:)= a_mtp_opt;
end
vA_opt = reshape(w_opt(starti:starti+NMuscle*N-1),NMuscle,N)';
starti = starti + NMuscle*N;
e_a_opt = reshape(w_opt(starti:starti+nq.arms*N-1),nq.arms,N)';
starti = starti + nq.arms*N;
if strcmp(S.subject.mtp_type,'active')
    e_mtp_opt = reshape(w_opt(starti:starti+nq.mtp*N-1),nq.mtp,N)';
    starti = starti + nq.mtp*N;
end
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
for k=1:N
    rangei = k*(d+1)-(d-1):k*(d+1);
    rangebi = (k-1)*d+1:k*d;
    a_mesh_col_opt(rangei,:) = a_col_opt(rangebi,:);
    FTtilde_mesh_col_opt(rangei,:) = FTtilde_col_opt(rangebi,:);
    Qs_mesh_col_opt(rangei,:) = Qs_col_opt(rangebi,:);
    Qdots_mesh_col_opt(rangei,:) = Qdots_col_opt(rangebi,:);
    a_a_mesh_col_opt(rangei,:) = a_a_col_opt(rangebi,:);
    if strcmp(S.subject.mtp_type,'active')
        a_mtp_mesh_col_opt(rangei,:) = a_mtp_col_opt(rangebi,:);
    end
end

%% Unscale results
% States at mesh points
% Qs (1:N-1)
q_opt_unsc.rad = Qs_opt(1:end-1,:).*repmat(...
    scaling.Qs,size(Qs_opt(1:end-1,:),1),1);
% Convert in degrees
q_opt_unsc.deg = q_opt_unsc.rad;
q_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations) = q_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations).*180/pi;
% Qs (1:N)
q_opt_unsc_all.rad = Qs_opt.*repmat(scaling.Qs,size(Qs_opt,1),1);
% Convert in degrees
q_opt_unsc_all.deg = q_opt_unsc_all.rad;
q_opt_unsc_all.deg(:,model_info.ExtFunIO.jointi.rotations) = q_opt_unsc_all.deg(:,model_info.ExtFunIO.jointi.rotations).*180/pi;
% Qdots (1:N-1)
qdot_opt_unsc.rad = Qdots_opt(1:end-1,:).*repmat(...
    scaling.Qdots,size(Qdots_opt(1:end-1,:),1),1);
% Convert in degrees
qdot_opt_unsc.deg = qdot_opt_unsc.rad;
qdot_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations) = qdot_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations).*180/pi;
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
if strcmp(S.subject.mtp_type,'active')
    % Mtp activations (1:N-1)
    a_mtp_opt_unsc = a_mtp_opt(1:end-1,:);
    % Mtp activations (1:N)
    a_mtp_opt_unsc_all = a_mtp_opt;
    % Arm activations
    a_a_col_opt_unsc = a_a_col_opt;
    % Mtp activations
    a_mtp_col_opt_unsc = a_mtp_col_opt;
    
    % Mtp excitations (control)
    e_mtp_opt_unsc = e_mtp_opt;
end
% Controls at mesh points
% Time derivative of muscle activations (states)
vA_opt_unsc = vA_opt.*repmat(scaling.vA,size(vA_opt,1),size(vA_opt,2));
% Get muscle excitations from time derivative of muscle activations
e_opt_unsc = computeExcitationRaasch(a_opt_unsc,vA_opt_unsc,...
    ones(1,NMuscle)*tdeact,ones(1,NMuscle)*tact);
% Arm excitations
e_a_opt_unsc = e_a_opt;
% States at collocation points
% Qs
q_col_opt_unsc.rad = Qs_col_opt.*repmat(scaling.Qs,size(Qs_col_opt,1),1);
% Convert in degrees
q_col_opt_unsc.deg = q_col_opt_unsc.rad;
q_col_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations) = q_col_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations).*180/pi;
% Qdots
qdot_col_opt_unsc.rad = Qdots_col_opt.*repmat(...
    scaling.Qdots,size(Qdots_col_opt,1),1);
% Convert in degrees
qdot_col_opt_unsc.deg = qdot_col_opt_unsc.rad;
qdot_col_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations) = qdot_col_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations).*180/pi;
% Muscle activations
a_col_opt_unsc = a_col_opt.*repmat(...
    scaling.a,size(a_col_opt,1),size(a_col_opt,2));
% Muscle-tendon forces
FTtilde_col_opt_unsc = FTtilde_col_opt.*repmat(...
    scaling.FTtilde,size(FTtilde_col_opt,1),1);
% "Slack" controls at collocation points
% Time derivative of Qdots
qdotdot_col_opt_unsc.rad = ...
    qdotdot_col_opt.*repmat(scaling.Qdotdots,size(qdotdot_col_opt,1),1);
% Convert in degrees
qdotdot_col_opt_unsc.deg = qdotdot_col_opt_unsc.rad;
qdotdot_col_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations) = qdotdot_col_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations).*180/pi;
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

%% Joint torques and ground reaction forces at mesh points (N), except #1
Xk_Qs_Qdots_opt             = zeros(N,2*nq.all);
Xk_Qs_Qdots_opt(:,1:2:end)  = q_opt_unsc_all.rad(2:end,:);
Xk_Qs_Qdots_opt(:,2:2:end)  = qdot_opt_unsc_all.rad(2:end,:);
Xk_Qdotdots_opt             = qdotdot_col_opt_unsc.rad(d:d:end,:);
Foutk_opt                   = zeros(N,F.nnz_out);
Tau_passk_opt_all           = zeros(N,nq.muscleActuated+nq.arms+nq.mtp);

for i = 1:N
    % ID moments
    [res] = F([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)']);
    Foutk_opt(i,:) = full(res);
    % passive moments
    Tau_passk_opt_all(i,:) = full(f_casadi.AllPassiveTorques(q_opt_unsc_all.rad(i+1,:),qdot_opt_unsc_all.rad(i+1,:)));
end
GRFk_opt = Foutk_opt(:,[model_info.ExtFunIO.GRFs.right_foot model_info.ExtFunIO.GRFs.left_foot]);

%% Joint torques and ground reaction forces at collocation points
Xj_Qs_Qdots_opt              = zeros(d*N,2*nq.all);
Xj_Qs_Qdots_opt(:,1:2:end)   = q_col_opt_unsc.rad;
Xj_Qs_Qdots_opt(:,2:2:end)   = qdot_col_opt_unsc.rad;
Xj_Qdotdots_opt              = qdotdot_col_opt_unsc.rad;
Foutj_opt                    = zeros(d*N,F.nnz_out);
Tau_passj_opt_all            = zeros(d*N,nq.muscleActuated+nq.arms+nq.mtp);
Tau_passj_opt_arm            = zeros(d*N,nq.arms);
Tau_passj_opt_mtp            = zeros(d*N,nq.mtp);
Tau_passj_opt_noMTP          = zeros(d*N,nq.muscleActuated+nq.arms);
Tau_passj_opt_muscleActuated = zeros(d*N,nq.muscleActuated);
for i = 1:d*N
    % inverse dynamics
    [res] = F([Xj_Qs_Qdots_opt(i,:)';Xj_Qdotdots_opt(i,:)']);
    Foutj_opt(i,:) = full(res);
    % passive torques
    Tau_passj_opt_all(i,:) = full(f_casadi.AllPassiveTorques(q_col_opt_unsc.rad(i,:),qdot_col_opt_unsc.rad(i,:)));
    Tau_passj_opt_arm(i,:) = Tau_passj_opt_all(i,(nq.muscleActuated+1):(nq.muscleActuated+nq.arms));
    Tau_passj_opt_mtp(i,:) = Tau_passj_opt_all(i,(nq.muscleActuated+nq.arms+1):end);
    Tau_passj_opt_noMTP(i,:) = Tau_passj_opt_all(i,1:(nq.muscleActuated+nq.arms));
    Tau_passj_opt_muscleActuated(i,:) = Tau_passj_opt_all(i,1:nq.muscleActuated);
end

%% Assert average speed
dist_trav_opt = q_opt_unsc_all.rad(end,model_info.ExtFunIO.coordi.pelvis_tx) - ...
    q_opt_unsc_all.rad(1,model_info.ExtFunIO.coordi.pelvis_tx); % distance traveled
time_elaps_opt = tf_opt; % time elapsed
vel_aver_opt = dist_trav_opt/time_elaps_opt;
% assert_v_tg should be 0
assert_v_tg = abs(vel_aver_opt-S.subject.v_pelvis_x_trgt);
if assert_v_tg > 1*10^(-S.solver.tol_ipopt)
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
vA_cost         = 0;
dFTtilde_cost   = 0;
QdotdotArm_cost = 0;
count           = 1;
h_opt           = tf_opt/N;
for k=1:N
    for j=1:d
        % Get muscle-tendon lengths, velocities, moment arms
        qin_opt_all = Xj_Qs_Qdots_opt(count,1:2:end);
        qdotin_opt_all = Xj_Qs_Qdots_opt(count,2:2:end);
        [lMTkj_opt_all,vMTkj_opt_all,~] = ...
            f_casadi.lMT_vMT_dM(qin_opt_all,qdotin_opt_all);
        % force equilibirum
        [~,~,Fce_opt_all,Fpass_opt_all,Fiso_opt_all] = ...
            f_casadi.forceEquilibrium_FtildeState_all_tendon(...
            a_col_opt_unsc(count,:)',FTtilde_col_opt_unsc(count,:)',...
            dFTtilde_col_opt_unsc(count,:)',full(lMTkj_opt_all),...
            full(vMTkj_opt_all),tensions);
        % muscle-tendon kinematics
        [~,lMtilde_opt_all] = f_casadi.FiberLength_TendonForce_tendon(...
            FTtilde_col_opt_unsc(count,:)',full(lMTkj_opt_all));
        [vM_opt_all,~] = f_casadi.FiberVelocity_TendonForce_tendon(...
            FTtilde_col_opt_unsc(count,:)',...
            dFTtilde_col_opt_unsc(count,:)',full(lMTkj_opt_all),...
            full(vMTkj_opt_all));
        
        % Bhargava et al. (2004)
        if strcmp(S.metabolicE.model,'Bhargava2004')
            [e_tot_all,~,~,~,~,~] = f_casadi.getMetabolicEnergySmooth2004all(...
            a_col_opt_unsc(count,:)',a_col_opt_unsc(count,:)',...
            full(lMtilde_opt_all),...
            full(vM_opt_all),full(Fce_opt_all),full(Fpass_opt_all),...
            MuscleMass',pctsts,full(Fiso_opt_all),S.subject.mass,10);
        else
            error('No energy model selected');
        end
        e_tot_opt_all = full(e_tot_all)';
        
        % objective function
        J_opt = J_opt + 1/(dist_trav_opt)*(...
            W.E*B(j+1)          *(f_casadi.J_N_muscles_exp(e_tot_opt_all,W.E_exp))/S.subject.mass*h_opt + ...
            W.a*B(j+1)          *(f_casadi.J_N_muscles(a_col_opt(count,:)))*h_opt + ...
            W.e_arm*B(j+1)      *(f_casadi.J_arms_dof(e_a_opt(k,:)))*h_opt +...
            W.q_dotdot*B(j+1)   *(f_casadi.J_noarms_dof(qdotdot_col_opt(count,model_info.ExtFunIO.jointi.noarmsi)))*h_opt + ...
            W.pass_torq*B(j+1)  *(f_casadi.J_pass_noMTP_dof(Tau_passj_opt_noMTP(count,:)))*h_opt + ... % MTP is not a part of this, Antoine had it in there
            W.slack_ctrl*B(j+1) *(f_casadi.J_N_muscles(vA_opt(k,:)))*h_opt + ...
            W.slack_ctrl*B(j+1) *(f_casadi.J_N_muscles(dFTtilde_col_opt(count,:)))*h_opt + ...
            W.slack_ctrl*B(j+1) *(f_casadi.J_arms_dof(qdotdot_col_opt(count,model_info.ExtFunIO.jointi.armsi)))*h_opt);

        if strcmp(S.subject.mtp_type,'active')
            J_opt = J_opt + 1/(dist_trav_opt)*(...
                W.e_mtp*B(j+1)  *(f_casadi.J_2(e_mtp_opt(k,:)))*h_opt);
            Mtp_cost = Mtp_cost + W.e_mtp*B(j+1)*...
                (f_casadi.J_2(e_mtp_opt(k,:)))*h_opt;
        end
        
        E_cost = E_cost + W.E*B(j+1)*...
            (f_casadi.J_N_muscles_exp(e_tot_opt_all,W.E_exp))/S.subject.mass*h_opt;
        A_cost = A_cost + W.a*B(j+1)*...
            (f_casadi.J_N_muscles(a_col_opt(count,:)))*h_opt;
        Arm_cost = Arm_cost + W.e_arm*B(j+1)*...
            (f_casadi.J_arms_dof(e_a_opt(k,:)))*h_opt;
        Qdotdot_cost = Qdotdot_cost + W.q_dotdot*B(j+1)*...
            (f_casadi.J_noarms_dof(qdotdot_col_opt(count,model_info.ExtFunIO.jointi.noarmsi)))*h_opt;
        Pass_cost = Pass_cost + W.pass_torq*B(j+1)*...
            (f_casadi.J_pass_noMTP_dof(Tau_passj_opt_noMTP(count,:)))*h_opt;
        vA_cost = vA_cost + W.slack_ctrl*B(j+1)*...
            (f_casadi.J_N_muscles(vA_opt(k,:)))*h_opt;
        dFTtilde_cost = dFTtilde_cost + W.slack_ctrl*B(j+1)*...
            (f_casadi.J_N_muscles(dFTtilde_col_opt(count,:)))*h_opt;
        QdotdotArm_cost = QdotdotArm_cost + W.slack_ctrl*B(j+1)*...
            (f_casadi.J_arms_dof(qdotdot_col_opt(count,model_info.ExtFunIO.jointi.armsi)))*h_opt;
        count = count + 1;
    end
end
J_optf = full(J_opt);           Obj.J = J_optf;
E_costf = full(E_cost);         Obj.E = E_costf;
A_costf = full(A_cost);         Obj.A = A_costf;
Arm_costf = full(Arm_cost);     Obj.Arm = Arm_costf;
Qdotdot_costf = full(Qdotdot_cost); Obj.qdd = Qdotdot_costf;
Pass_costf = full(Pass_cost);   Obj.Pass = Pass_costf;
vA_costf = full(vA_cost);       Obj.vA = vA_costf;
dFTtilde_costf = full(dFTtilde_cost); Obj.dFTtilde = dFTtilde_costf;
QdotdotArm_costf = full(QdotdotArm_cost); Obj.qdd_arm = QdotdotArm_costf;
if strcmp(S.subject.mtp_type,'active')
    Mtp_costf = full(Mtp_cost);     Obj.Mtp = Mtp_costf;
    contributionCost.absoluteValues = 1/(dist_trav_opt)*[E_costf,A_costf,...
        Arm_costf,Qdotdot_costf,Pass_costf,vA_costf,dFTtilde_costf,...
        QdotdotArm_costf,Mtp_costf];
    contributionCost.relativeValues = 1/(dist_trav_opt)*[E_costf,A_costf,...
        Arm_costf,Qdotdot_costf,Pass_costf,vA_costf,dFTtilde_costf,...
        QdotdotArm_costf,Mtp_costf]./J_optf*100;
    contributionCost.relativeValuesRound2 = ...
        round(contributionCost.relativeValues,2);
    contributionCost.labels = {'metabolicEnergy','muscleActivation',...
        'armExcitation','jointAccelerations','passiveTorques','dadt','dFdt',...
        'armAccelerations','mtpExcitation'};
    % assertCost should be 0
    assertCost = abs(J_optf - 1/(dist_trav_opt)*(E_costf+A_costf+Arm_costf+...
        Mtp_costf+Qdotdot_costf+Pass_costf+vA_costf+dFTtilde_costf+...
        QdotdotArm_costf));
else
    contributionCost.absoluteValues = 1/(dist_trav_opt)*[E_costf,A_costf,...
        Arm_costf,Qdotdot_costf,Pass_costf,vA_costf,dFTtilde_costf,...
        QdotdotArm_costf];
    contributionCost.relativeValues = 1/(dist_trav_opt)*[E_costf,A_costf,...
        Arm_costf,Qdotdot_costf,Pass_costf,vA_costf,dFTtilde_costf,...
        QdotdotArm_costf]./J_optf*100;
    contributionCost.relativeValuesRound2 = ...
        round(contributionCost.relativeValues,2);
    contributionCost.labels = {'metabolicEnergy','muscleActivation',...
        'armExcitation','jointAccelerations','passiveTorques','dadt','dFdt',...
        'armAccelerations'};
    % assertCost should be 0
    assertCost = abs(J_optf - 1/(dist_trav_opt)*(E_costf+A_costf+Arm_costf+...
        Qdotdot_costf+Pass_costf+vA_costf+dFTtilde_costf+...
        QdotdotArm_costf));
end    
assertCost2 = abs(stats.iterations.obj(end) - J_optf);
if assertCost > 1*10^(-S.solver.tol_ipopt)
    disp('Issue when reconstructing optimal cost wrt sum of terms')
end
if assertCost2 > 1*10^(-S.solver.tol_ipopt)
    disp('Issue when reconstructing optimal cost wrt stats')
end
Outname = fullfile(OutFolder,[S.subject.name '_ContributionCost.mat']);
save(Outname,'contributionCost');

%% Reconstruct full gait cycle
if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    % detect heelstrike
    [IC1i_c,IC1i_s,HS1] = getHeelstrikeSimulation(GRFk_opt,N);
        
    % Qs
    Qs_GC = zeros(N*2,size(q_opt_unsc.deg,2));
    Qs_GC(1:N-IC1i_s+1,:) = q_opt_unsc.deg(IC1i_s:end,:);
    Qs_GC(N-IC1i_s+2:N-IC1i_s+1+N,QsSymA) = q_opt_unsc.deg(1:end,QsSymB);
    Qs_GC(N-IC1i_s+2:N-IC1i_s+1+N,QsOpp) = -q_opt_unsc.deg(1:end,QsOpp);
    Qs_GC(N-IC1i_s+2:N-IC1i_s+1+N,model_info.ExtFunIO.coordi.pelvis_tx) = ...
        q_opt_unsc.deg(1:end,model_info.ExtFunIO.coordi.pelvis_tx) + ...
        q_opt_unsc_all.deg(end,model_info.ExtFunIO.coordi.pelvis_tx);
    Qs_GC(N-IC1i_s+2+N:2*N,:) = q_opt_unsc.deg(1:IC1i_s-1,:);
    Qs_GC(N-IC1i_s+2+N:2*N,model_info.ExtFunIO.coordi.pelvis_tx) = ...
        q_opt_unsc.deg(1:IC1i_s-1,model_info.ExtFunIO.coordi.pelvis_tx) + ...
        2*q_opt_unsc_all.deg(end,model_info.ExtFunIO.coordi.pelvis_tx);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        Qs_GC(:,QsSymA_ptx)  = Qs_GC(:,QsSymB_ptx);
        Qs_GC(:,QsOpp)       = -Qs_GC(:,QsOpp);
    end
    temp_Qs_GC_pelvis_tx = Qs_GC(1,model_info.ExtFunIO.coordi.pelvis_tx);
    Qs_GC(:,model_info.ExtFunIO.coordi.pelvis_tx) = Qs_GC(:,model_info.ExtFunIO.coordi.pelvis_tx)-...
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
    Ts_opt = Ts_opt./S.subject.mass;
    
    % Muscle-Tendon Forces
    orderMusInv = model_info.ExtFunIO.symQs.orderMusInv;
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
    orderArmInv = model_info.ExtFunIO.symQs.orderArmInv;
    a_a_GC = zeros(N*2,nq.arms);
    a_a_GC(1:N-IC1i_s+1,:) = a_a_opt_unsc(IC1i_s:end,:);
    a_a_GC(N-IC1i_s+2:N-IC1i_s+1+N,:) = a_a_opt_unsc(1:end,orderArmInv);
    a_a_GC(N-IC1i_s+2+N:2*N,:) = a_a_opt_unsc(1:IC1i_s-1,:);
    % If the first heel strike was on the left foot then we invert so that
    % we always start with the right foot, for analysis purpose
    if strcmp(HS1,'l')
        a_a_GC(:,:) = a_a_GC(:,orderArmInv);
    end
    
    if strcmp(S.subject.mtp_type,'active')
        % Mtp activations
        orderMtpInv = [2 1]; % There are only 2 mtps
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
    end    
    % Passive joint torques
    Tau_pass_opt_inv = model_info.ExtFunIO.symQs.orderTauPassInv;
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
    q_opt_GUI_GC = zeros(2*N,1+nq.all+2); % added pro_sup for both arms
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

JointAngle.labels = [{'time'}, fields(model_info.ExtFunIO.coordi)', {'pro_sup_l'},{'pro_sup_r'}];
% Two gait cycles
% Joint angles
q_opt_GUI_GC_2 = [q_opt_GUI_GC;q_opt_GUI_GC];
q_opt_GUI_GC_2(2*N+1:4*N,1) = q_opt_GUI_GC_2(2*N+1:4*N,1) + ...
    q_opt_GUI_GC_2(end,1) + ...
    q_opt_GUI_GC_2(end,1)-q_opt_GUI_GC_2(end-1,1);
q_opt_GUI_GC_2(2*N+1:4*N,model_info.ExtFunIO.coordi.pelvis_tx+1) = ...
    q_opt_GUI_GC_2(2*N+1:4*N,model_info.ExtFunIO.coordi.pelvis_tx+1) + ...
    2*q_opt_unsc_all.deg(end,model_info.ExtFunIO.coordi.pelvis_tx);
% Muscle activations (to have muscles turning red when activated).
Acts_GC_GUI = [Acts_GC;Acts_GC];
% Combine data joint angles and muscle activations
JointAngleMuscleAct.data = [q_opt_GUI_GC_2,Acts_GC_GUI];
% Get muscle labels
muscleNamesAll = model_info.muscle_info.muscle_names;
% Combine labels joint angles and muscle activations
JointAngleMuscleAct.labels = JointAngle.labels;
for i = 1:NMuscle
    JointAngleMuscleAct.labels{i+size(q_opt_GUI_GC_2,2)} = ...
        [muscleNamesAll{i},'/activation'];
end
JointAngleMuscleAct.inDeg = 'yes';
filenameJointAngles = fullfile(OutFolder,[S.subject.name '.mot']);
write_motionFile_v40(JointAngleMuscleAct, filenameJointAngles);

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
R.Tid       = Ts_opt.*S.subject.mass;
R.a         = Acts_GC;
R.e         = e_GC;
R.S           = S;  % settings for post processing
R.Sopt        = Sopt; % original settings used to solve the OCP
R.body_mass   = S.subject.mass;
R.a_arm       = a_a_GC;
R.e_arm       = e_a_GC;
R.a_mtp       = a_mtp_GC;
R.e_mtp       = e_mtp_GC;
R.TPass       = Tau_pass_opt_GC;
R.Obj         = Obj;
R.BodyKin     = BodyKin;

FilenameAnalysis = fullfile(OutFolder,[S.subject.name '_pp.mat']);
save(FilenameAnalysis,'R');

end

