function [] = OCP_formulation(S,model_info,f_casadi)
% --------------------------------------------------------------------------
% OCP_formulation
%   This function formulates the OCP and calls the solver
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
%   - f_casadi -
%   * Struct containing all casadi functions.
%
% OUTPUT:
%   - This function returns no outputs -
% 
% Original author: Dhruv Gupta and Lars D'Hondt
% Original date: January-May/2022
%
% Last edit by: Bram Van Den Bosch
% Last edit date: 23/Sept/2024
% --------------------------------------------------------------------------

disp('Start formulating OCP...')
t0 = tic;

%% User inputs (typical settings structure)
% settings for optimization
N = S.solver.N_meshes; % number of mesh intervals
W = S.weights; % weights optimization
nq = model_info.ExtFunIO.jointi.nq; % lengths of coordinate subsets

%% Load external functions
import casadi.*
% The external function performs inverse dynamics through the
% OpenSim/Simbody C++ API. This external function is compiled as a dll from
% which we create a Function instance using CasADi in MATLAB. More details
% about the external function can be found in the documentation.
pathmain = pwd;
% [filepath,~,~] = fileparts(mfilename('fullpath'));
% [pathRepo,~,~] = fileparts(filepath);
% addpath(genpath(pathRepo));
% Loading external functions.
setup.derivatives =  'AD'; % Algorithmic differentiation
cd(S.misc.subject_path)
F  = external('F', fullfile(S.misc.subject_path, S.misc.external_function));
cd(pathmain);

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
    model_info.muscle_info.muscle_spanning_joint_info);
% calculate total number of joints that each muscle crosses (used later)
sumCross = sum(model_info.muscle_info.muscle_spanning_joint_info);

% Parameters for activation dynamics
tact = model_info.muscle_info.tact; % Activation time constant
tdeact = model_info.muscle_info.tdeact; % Deactivation time constant

% Muscles indices for reading synergies
if (S.subject.synergies) 
    idx_m_r = model_info.muscle_info.idx_right;
    idx_m_l = model_info.muscle_info.idx_left;
    muscleNames_r = muscleNames(idx_m_r);
    muscleNames_l = muscleNames(idx_m_l);
end

%% Metabolic energy model parameters
% We extract the specific tensions and slow twitch rations.
tensions = struct_array_to_double_array(model_info.muscle_info.parameters,'specific_tension');
pctsts = struct_array_to_double_array(model_info.muscle_info.parameters,'slow_twitch_fiber_ratio');

%% Function to compute muscle mass
MuscleMass = struct_array_to_double_array(model_info.muscle_info.parameters,'muscle_mass');

%% Get bounds and initial guess

bounds_nsc = getBounds(S,model_info);
scaling = getScaleFactor(S,model_info,bounds_nsc);
bounds = scaleBounds(S,model_info,bounds_nsc,scaling);

if strcmp(S.solver.IG_selection,'quasi-random')
    guess = getGuess_QR_opti(S,model_info,scaling,d);
else
    guess = getGuess_DI_opti(S,model_info,scaling,d);
end


if (S.misc.visualize_bounds)
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
lboundsQsk(model_info.ExtFunIO.jointi.base_forward,1) = ...
    bounds.Qs_0.lower(model_info.ExtFunIO.jointi.base_forward);
uboundsQsk = bounds.Qs.upper'*ones(1,N+1);
uboundsQsk(model_info.ExtFunIO.jointi.base_forward,1) = ...
    bounds.Qs_0.upper(model_info.ExtFunIO.jointi.base_forward);
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
% Torque actuator activations at mesh points
if nq.torqAct > 0
    a_a = opti.variable(nq.torqAct,N+1);
    opti.subject_to(bounds.a_a.lower'*ones(1,N+1) < a_a < ...
        bounds.a_a.upper'*ones(1,N+1));
    opti.set_initial(a_a, guess.a_a');
    % Torque actuator activations at collocation points
    a_a_col = opti.variable(nq.torqAct,d*N);
    opti.subject_to(bounds.a_a.lower'*ones(1,d*N) < a_a_col < ...
        bounds.a_a.upper'*ones(1,d*N));
    opti.set_initial(a_a_col, guess.a_a_col');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define controls
% Time derivative of muscle activations (states) at mesh points
vA = opti.variable(NMuscle, N);
opti.subject_to(bounds.vA.lower'*ones(1,N) < vA < ...
    bounds.vA.upper'*ones(1,N));
opti.set_initial(vA, guess.vA');
% Torque actuator excitations
if nq.torqAct > 0
    e_a = opti.variable(nq.torqAct, N);
    opti.subject_to(bounds.e_a.lower'*ones(1,N) < e_a < ...
        bounds.e_a.upper'*ones(1,N));
    opti.set_initial(e_a, guess.e_a');
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

%% Helper function for orthoses
% variables
a_MX = MX.sym('a',NMuscle,N);
Qs_MX = MX.sym('Qs',nq.all,N);
Qdots_MX = MX.sym('Qdots',nq.all,N);
Qddots_MX = MX.sym('Qddots',nq.all,N);

% unscale variables
Qs_MX_nsc = Qs_MX.*(scaling.Qs'*ones(1,size(Qs_MX,2)));
Qdots_MX_nsc = Qdots_MX.*(scaling.Qdots'*ones(1,size(Qdots_MX,2)));
Qddots_MX_nsc = Qddots_MX.*(scaling.Qdotdots'*ones(1,size(Qddots_MX,2)));

% evaluate orthosis function
[M_ort_coord_MX, M_ort_body_MX] = f_casadi.f_orthosis_mesh_all(Qs_MX_nsc, Qdots_MX_nsc,...
    Qddots_MX_nsc, a_MX);

% create function
f_orthosis_mesh_all = Function('f_orthosis_mesh_all',{Qs_MX, Qdots_MX,...
    Qddots_MX, a_MX},{M_ort_coord_MX, M_ort_body_MX});

% evaluate helper function
[M_ort_coord_opti, M_ort_body_opti] = f_orthosis_mesh_all(Qs(:,1:N), Qdots(:,1:N),...
    A_col(:,1:3:3*N), a(:,1:N)); 
    % note: A is at 1st collocation point of mesh interval instead of at 1st mesh point


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If Muscle synergies, add additional variables:
    % Synergy activations as states at mesh points (right and left independently)
    % Synergy weights as static parameters (constant in time)
if (S.subject.synergies)        
    if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle') % Same number of synergies right and left, 
        % and same weights right and left

        % Right synergy activations at mesh points
        SynH_r = opti.variable(S.subject.NSyn_r,N+1);
        opti.subject_to(bounds.SynH.lower(1:S.subject.NSyn_r)'*ones(1,N+1) < SynH_r < ...
            bounds.SynH.upper(1:S.subject.NSyn_r)'*ones(1,N+1));
        opti.set_initial(SynH_r, guess.SynH(:,1:S.subject.NSyn_r)');

        % Left synergy activations at mesh points
        SynH_l = opti.variable(S.subject.NSyn_l,N+1);
        opti.subject_to(bounds.SynH.lower(1:S.subject.NSyn_l)'*ones(1,N+1) < SynH_l < ...
            bounds.SynH.upper(1:S.subject.NSyn_l)'*ones(1,N+1));
        opti.set_initial(SynH_l, guess.SynH(:,1:S.subject.NSyn_l)');

        % Synergy weights
        SynW_r = opti.variable(length(idx_m_r),S.subject.NSyn_r);
        opti.subject_to(bounds.SynW.lower*ones(length(idx_m_r),S.subject.NSyn_r) < SynW_r < bounds.SynW.upper*ones(length(idx_m_r),S.subject.NSyn_r));
        opti.set_initial(SynW_r,  guess.SynW*ones(length(idx_m_r),S.subject.NSyn_r));

        SynW_l = SynW_r; 

    elseif strcmp(S.misc.gaitmotion_type,'FullGaitCycle')

        % Right synergy activations at mesh points
        SynH_r = opti.variable(S.subject.NSyn_r,N+1);
        % Bounds and initial guess from muscle activations (this could be
        % modified)
        opti.subject_to(bounds.SynH.lower(1:S.subject.NSyn_r)'*ones(1,N+1) < SynH_r < ...
            bounds.SynH.upper(1:S.subject.NSyn_r)'*ones(1,N+1));
        opti.set_initial(SynH_r, guess.SynH(:,1:S.subject.NSyn_r)');

        % Left synergy activations at mesh points (left)
        SynH_l = opti.variable(S.subject.NSyn_l,N+1);
        opti.subject_to(bounds.SynH.lower(1:S.subject.NSyn_l)'*ones(1,N+1) < SynH_l < ...
            bounds.SynH.upper(1:S.subject.NSyn_l)'*ones(1,N+1));
        opti.set_initial(SynH_l, guess.SynH(:,1:S.subject.NSyn_l)');

        % Right synergy weights
        SynW_r = opti.variable(length(idx_m_r),S.subject.NSyn_r);
        opti.subject_to(bounds.SynW.lower*ones(length(idx_m_r),S.subject.NSyn_r) < SynW_r < bounds.SynW.upper*ones(length(idx_m_r),S.subject.NSyn_r));
        opti.set_initial(SynW_r, guess.SynW*ones(length(idx_m_r),S.subject.NSyn_r));
        
        % Left synergy weights
        SynW_l = opti.variable(length(idx_m_l),S.subject.NSyn_l);
        opti.subject_to(bounds.SynW.lower*ones(length(idx_m_l),S.subject.NSyn_l) < SynW_l < bounds.SynW.upper*ones(length(idx_m_l),S.subject.NSyn_l));
        opti.set_initial(SynW_l, guess.SynW*ones(length(idx_m_l),S.subject.NSyn_l));
    end
end

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
if nq.torqAct > 0
    a_ak        = MX.sym('a_ak',nq.torqAct);
    a_aj        = MX.sym('a_akmesh',nq.torqAct,d);
    a_akj       = [a_ak a_aj];
end
% Define CasADi variables for controls
vAk     = MX.sym('vAk',NMuscle);
if nq.torqAct > 0
    e_ak    = MX.sym('e_ak',nq.torqAct);
end

% Define CasADi variables for "slack" controls
dFTtildej   = MX.sym('dFTtildej',NMuscle,d);
Aj          = MX.sym('Aj',nq.all,d);

% Define CasADi variables for orthosis moments (or forces)
M_ort_coordk = MX.sym('M_ort_coord',nq.all,1); % moments on coordinates
M_ort_bodyk = MX.sym('M_ort_body',model_info.ExtFunIO.input.nInputs,1); % moments on bodies

% If muscle synergies, define CasADi variables for additional variables
if (S.subject.synergies)
    SynH_rk         = MX.sym('SynH_rk',S.subject.NSyn_r);
    SynH_lk         = MX.sym('SynH_lk',S.subject.NSyn_l);
    SynW_rk         = MX.sym('SynW_rk',length(idx_m_r),S.subject.NSyn_r);
    SynW_lk         = MX.sym('SynW_lk',length(idx_m_l),S.subject.NSyn_l);
end

J           = 0; % Initialize cost function
eq_constr   = {}; % Initialize equality constraint vector
ineq_constr_deact = {}; % Initialize inequality constraint vector
ineq_constr_act = {}; % Initialize inequality constraint vector
ineq_constr_distance = cell(length(S.bounds.distanceConstraints),1);
ineq_constr_syn = {}; % Initialize inequality constraint vector

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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get muscle-tendon lengths, velocities, and moment arms
    qinj    = Qskj_nsc(:, j+1);
    qdotinj = Qdotskj_nsc(:, j+1);
    [lMTj,vMTj,MAj] =  f_casadi.lMT_vMT_dM(qinj',qdotinj');
    % Derive the moment arms of all the muscles crossing each joint
    for i=1:nq.musAct
        MA_j.(model_info.ExtFunIO.coord_names.muscleActuated{i}) = ...
            MAj(mai(model_info.ExtFunIO.jointi.muscleActuated(i)).mus',...
            model_info.ExtFunIO.jointi.muscleActuated(i));
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
            MuscleMass',pctsts,Fisoj,model_info.mass,S.metabolicE.tanh_b);
    else
        error('No energy model selected');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get passive joint torques for dynamics
    Tau_passj = f_casadi.AllPassiveTorques(Qskj_nsc(:,j+1),Qdotskj_nsc(:,j+1));
    % Get passive joint torques for cost function
    Tau_passj_cost = f_casadi.AllPassiveTorques_cost(Qskj_nsc(:,j+1),Qdotskj_nsc(:,j+1));
    
    % Expression for the state derivatives at the collocation points
    Qsp_nsc      = Qskj_nsc*C(:,j+1);
    Qdotsp_nsc   = Qdotskj_nsc*C(:,j+1);
    FTtildep_nsc = FTtildekj_nsc*C(:,j+1);
    ap           = akj*C(:,j+1);
    if nq.torqAct > 0
        a_ap         = a_akj*C(:,j+1);
    end
    % Append collocation equations
    % Dynamic constraints are scaled using the same scale
    % factors as the ones used to scale the states
    % Activation dynamics (implicit formulation)
    eq_constr{end+1} = (h*vAk_nsc - ap)./scaling.a;
    % Contraction dynamics (implicit formulation)
    eq_constr{end+1} = (h*dFTtildej_nsc(:,j) - FTtildep_nsc)./scaling.FTtilde';
    % Skeleton dynamics (implicit formulation)
    qdotj_nsc = Qdotskj_nsc(:,j+1); % velocity
    eq_constr{end+1} = (h*qdotj_nsc - Qsp_nsc)./scaling.Qs';
    eq_constr{end+1} = (h*Aj_nsc(:,j) - Qdotsp_nsc)./scaling.Qdots';
    % Torque actuator activation dynamics (explicit formulation)
    if nq.torqAct > 0
        da_adtj = f_casadi.ActuatorActivationDynamics(e_ak,a_akj(:,j+1));
        eq_constr{end+1} = (h*da_adtj - a_ap)./scaling.a_a;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Orthosis moments on collocation point
    [M_ort_coordj, M_ort_bodyj] = f_casadi.f_orthosis_mesh_k(Qskj_nsc(:,j+1),...
        Qdotskj_nsc(:,j+1), Aj_nsc(:,j), akj(:,j+1));

    % add orthosis moments from input variables
    M_ort_coord_totj = M_ort_coordk + M_ort_coordj;
    M_ort_body_totj = M_ort_bodyk + M_ort_bodyj;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Add contribution to the cost function
    J = J + ...
        W.E          * B(j+1) *(f_casadi.J_muscles_exp(e_totj, W.E_exp))/model_info.mass*h + ...
        W.a          * B(j+1) *(f_casadi.J_muscles_exp(akj(:,j+1)', W.a_exp))*h + ...
        W.q_dotdot   * B(j+1) *(f_casadi.J_not_arms_dof(Aj(model_info.ExtFunIO.jointi.noarmsi,j)))*h + ...
        W.pass_torq  * B(j+1) *(f_casadi.J_lim_torq(Tau_passj_cost))*h + ...
        W.slack_ctrl * B(j+1) *(f_casadi.J_muscles(vAk))*h + ...
        W.slack_ctrl * B(j+1) *(f_casadi.J_muscles(dFTtildej(:,j)))*h;

    if nq.torqAct > 0
        J = J + W.e_torqAct      * B(j+1) *(f_casadi.J_torq_act(e_ak))*h;
    end
    if nq.arms > 0
        J = J + W.slack_ctrl * B(j+1) *(f_casadi.J_arms_dof(Aj(model_info.ExtFunIO.jointi.armsi,j)))*h;
    end

    % If muscle synergies: Instead of a - WH = 0 as an equality constraint, have it as a
        % term in the cost function to be minimized (+ inequality constraint)     
    if (S.subject.synergies)        
            syn_constr_k_r = ak(idx_m_r) - SynW_rk*SynH_rk;
            syn_constr_k_l = ak(idx_m_l) - SynW_lk*SynH_lk;
        J = J + W.SynConstr * B(j+1) *(f_casadi.J_muscles([syn_constr_k_r;syn_constr_k_l]))*h;   
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create zero (sparse) input vector for external function
    F_ext_input = MX(model_info.ExtFunIO.input.nInputs,1);
    % Assign Qs
    F_ext_input(model_info.ExtFunIO.input.Qs.all,1) = Qskj_nsc(:,j+1);
    % Assign Qdots
    F_ext_input(model_info.ExtFunIO.input.Qdots.all,1) = Qdotskj_nsc(:,j+1);
    % Assign Qdotdots (A)
    F_ext_input(model_info.ExtFunIO.input.Qdotdots.all,1) = Aj_nsc(:,j);
    % Assign forces and moments 
    F_ext_input = F_ext_input + M_ort_body_totj; % body forces and body moments from orthoses

    % Evaluate external function
    [Tj] = F(F_ext_input);

    % Evaluate ligament moment
    M_lig_j = f_casadi.ligamentMoment(Qskj_nsc(:,j+1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add path constraints
    for i=1:nq.all
        % total coordinate torque
        Ti = 0;

        % muscle moment
        cond_special_mtp = strcmp(S.subject.mtp_type,'2022paper') &&...
            contains(model_info.ExtFunIO.coord_names.all{i},'mtp');

        if ismember(i,model_info.ExtFunIO.jointi.muscleActuated) && ~cond_special_mtp
            % muscle forces
            FTj_coord_i = FTj(mai(i).mus',1);
            % total muscle moment
            M_mus_i = f_casadi.(['musc_cross_' num2str(sumCross(i))])...
            (MA_j.(model_info.ExtFunIO.coord_names.all{i}),FTj_coord_i);
            % add to total moment
            Ti = Ti + M_mus_i;
        end
        
        % torque actuator
        if nq.torqAct > 0 && ismember(i,model_info.ExtFunIO.jointi.torqueActuated)
            idx_act_i = find(model_info.ExtFunIO.jointi.torqueActuated(:)==i);
            T_act_i = a_akj(idx_act_i,j+1).*scaling.ActuatorTorque(idx_act_i);
            Ti = Ti + T_act_i;
        end

        % ligament moment
        Ti = Ti + M_lig_j(i);
        
        % passive moment
        if ~ismember(i,model_info.ExtFunIO.jointi.floating_base)
            Ti = Ti + Tau_passj(i);
        end
        
        % orthosis
        Ti = Ti + M_ort_coord_totj(i);

        % total coordinate torque equals inverse dynamics torque
        eq_constr{end+1} = (Tj(i,1) - Ti)./scaling.Moments(i);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Activation dynamics (implicit formulation)
    act1 = vAk_nsc + akj(:,j+1)./(ones(size(akj(:,j+1),1),1)*tdeact);
    act2 = vAk_nsc + akj(:,j+1)./(ones(size(akj(:,j+1),1),1)*tact);
    ineq_constr_deact{end+1} = act1;
    ineq_constr_act{end+1} = act2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contraction dynamics (implicit formulation)
    eq_constr{end+1} = Hilldiffj;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constraints to prevent parts of the skeleton to penetrate each other.
    for i_dc=1:length(S.bounds.distanceConstraints)
        % get position of points
        point1_i = S.bounds.points(strcmp({S.bounds.points.name},S.bounds.distanceConstraints(i_dc).point1));
        if strcmpi(point1_i.body,'ground')
            % position of point in ground is known
            pos1j = point1_i.point_in_body;
        else
            % position of point in body is provided by external function
            pos1j = Tj(model_info.ExtFunIO.position.(S.bounds.distanceConstraints(i_dc).point1),1);
        end
        point2_i = S.bounds.points(strcmp({S.bounds.points.name},S.bounds.distanceConstraints(i_dc).point2));
        if strcmpi(point2_i.body,'ground')
            % position of point in ground is known
            pos2j = point2_i.point_in_body;
        else
            % position of point in body is provided by external function
            pos2j = Tj(model_info.ExtFunIO.position.(S.bounds.distanceConstraints(i_dc).point2),1);
        end

        % calculate distance
        if length(S.bounds.distanceConstraints(i_dc).directionVectorIdx)==3 % 3D
            Qconstr = f_casadi.J_nn_3(pos1j(S.bounds.distanceConstraints(i_dc).directionVectorIdx) ...
                - pos2j(S.bounds.distanceConstraints(i_dc).directionVectorIdx));

            ineq_constr_distance{i_dc}{end+1} = Qconstr;

        elseif length(S.bounds.distanceConstraints(i_dc).directionVectorIdx)==2 % 2D
            Qconstr = f_casadi.J_nn_2(pos1j(S.bounds.distanceConstraints(i_dc).directionVectorIdx) ...
                - pos2j(S.bounds.distanceConstraints(i_dc).directionVectorIdx));

            ineq_constr_distance{i_dc}{end+1} = Qconstr;

        elseif length(S.bounds.distanceConstraints(i_dc).directionVectorIdx)==1 % 1D
            Qconstr = pos1j(S.bounds.distanceConstraints(i_dc).directionVectorIdx) ...
                - pos2j(S.bounds.distanceConstraints(i_dc).directionVectorIdx);

            ineq_constr_distance{i_dc}{end+1} = Qconstr;
        end

    end

end % End loop over collocation points

% Add tracking terms in the cost function if synergy weights are tracked
% Here we select the weights that we want to impose/track 
% (there are no conditions/constraints applied to the other weights)
if (S.subject.synergies) && (S.subject.TrackSynW)
    J_TrackSynW = W.TrackSynW*f_casadi.TrackSynW(SynW_rk, SynW_lk);
    J = J + J_TrackSynW;
end

% Synergies: a - WH = 0
% Only applied for mesh points
if (S.subject.synergies)
    ineq_constr_syn{end+1} =  [ak(idx_m_r);ak(idx_m_l)] - [SynW_rk*SynH_rk;SynW_lk*SynH_lk];
end

eq_constr = vertcat(eq_constr{:});
ineq_constr_deact = vertcat(ineq_constr_deact{:});
ineq_constr_act = vertcat(ineq_constr_act{:});
for i_dc=1:length(ineq_constr_distance)
    ineq_constr_distance{i_dc} = vertcat(ineq_constr_distance{i_dc}{:});
end
ineq_constr_syn = vertcat(ineq_constr_syn{:});

% Casadi function to get constraints and objective
coll_input_vars_def = {tfk,ak,aj,FTtildek,FTtildej,Qsk,Qsj,Qdotsk,Qdotsj,vAk,dFTtildej,Aj,M_ort_coordk,M_ort_bodyk};
if nq.torqAct > 0
    coll_input_vars_def = [coll_input_vars_def,{a_ak,a_aj,e_ak}];
end
if (S.subject.synergies)
    coll_input_vars_def = [coll_input_vars_def,{SynH_rk,SynH_lk,SynW_rk,SynW_lk}];
end

f_coll = Function('f_coll',coll_input_vars_def,...
        {eq_constr, ineq_constr_deact, ineq_constr_act,...
        ineq_constr_distance{:}, ineq_constr_syn,J});

% Repeat function for each mesh interval and assign evaluation to multiple threads
f_coll_map = f_coll.map(N,S.solver.parallel_mode,S.solver.N_threads);

% evaluate function with opti variables
coll_input_vars_eval = {tf,a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col,...
    Qs(:,1:end-1), Qs_col, Qdots(:,1:end-1), Qdots_col, vA, dFTtilde_col, A_col,...
     M_ort_coord_opti, M_ort_body_opti};
if nq.torqAct > 0
    coll_input_vars_eval = [coll_input_vars_eval, {a_a(:,1:end-1), a_a_col, e_a}];
end
if (S.subject.synergies)    
    coll_input_vars_eval = [coll_input_vars_eval,{SynH_r(:,1:end-1), SynH_l(:,1:end-1), SynW_r,SynW_l}];
end

coll_ineq_constr_distance = cell(1,length(ineq_constr_distance));

[coll_eq_constr, coll_ineq_constr_deact, coll_ineq_constr_act,...
    coll_ineq_constr_distance{:}, coll_ineq_constr_syn,Jall] =...
    f_coll_map(coll_input_vars_eval{:});

% equality constraints
opti.subject_to(coll_eq_constr == 0);

% inequality constraints (logical indexing not possible in MX arrays)
opti.subject_to(coll_ineq_constr_deact(:) >= 0);
opti.subject_to(coll_ineq_constr_act(:) <= 1/tact);

for i_dc=1:length(ineq_constr_distance)
    coll_ineq_constr_distance_i_dc = coll_ineq_constr_distance{i_dc};

    if ~isempty(S.bounds.distanceConstraints(i_dc).lower_bound) && ...
             ~isempty(S.bounds.distanceConstraints(i_dc).upper_bound)
        % lower and upper bound
        opti.subject_to(S.bounds.distanceConstraints(i_dc).lower_bound < coll_ineq_constr_distance_i_dc(:)...
            < S.bounds.distanceConstraints(i_dc).upper_bound);
    elseif ~isempty(S.bounds.distanceConstraints(i_dc).lower_bound)
        % only lower bound
        opti.subject_to(S.bounds.distanceConstraints(i_dc).lower_bound < coll_ineq_constr_distance_i_dc(:))
    elseif ~isempty(S.bounds.distanceConstraints(i_dc).upper_bound)
        % only upper bound
        opti.subject_to(coll_ineq_constr_distance_i_dc(:) < S.bounds.distanceConstraints(i_dc).upper_bound);
    end

end
if (S.subject.synergies)
    opti.subject_to(S.bounds.SynConstr.lower < coll_ineq_constr_syn(:) < S.bounds.SynConstr.upper); 
end

% Loop over mesh points
for k=1:N
    % Variables within current mesh interval
    % States
    akj = [a(:,k), a_col(:,(k-1)*d+1:k*d)];
    FTtildekj = [FTtilde(:,k), FTtilde_col(:,(k-1)*d+1:k*d)];
    Qskj = [Qs(:,k), Qs_col(:,(k-1)*d+1:k*d)];
    Qdotskj = [Qdots(:,k), Qdots_col(:,(k-1)*d+1:k*d)];
    if nq.torqAct > 0
        a_akj = [a_a(:,k), a_a_col(:,(k-1)*d+1:k*d)];
    end

    % Add equality constraints (next interval starts with end values of
    % states from previous interval)
    opti.subject_to(a(:,k+1) == akj*D);
    opti.subject_to(FTtilde(:,k+1) == FTtildekj*D); % scaled
    opti.subject_to(Qs(:,k+1) == Qskj*D); % scaled
    opti.subject_to(Qdots(:,k+1) == Qdotskj*D); % scaled
    if nq.torqAct > 0
        opti.subject_to(a_a(:,k+1) == a_akj*D);
    end

end % End loop over mesh points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional path constraints
if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    % Periodicity of the states (or rather LR symmetry -half gait cycle)
    % Qs and Qdots
    opti.subject_to(Qs(model_info.ExtFunIO.symQs.QsInvA,end) - Qs(model_info.ExtFunIO.symQs.QsInvB,1) == 0);
    opti.subject_to(Qdots(model_info.ExtFunIO.symQs.QdotsInvA,end) - Qdots(model_info.ExtFunIO.symQs.QdotsInvB,1) == 0);
    if ~isempty(model_info.ExtFunIO.symQs.QsOpp)
        opti.subject_to(Qs(model_info.ExtFunIO.symQs.QsOpp,end) + Qs(model_info.ExtFunIO.symQs.QsOpp,1) == 0);
        opti.subject_to(Qdots(model_info.ExtFunIO.symQs.QsOpp,end) + Qdots(model_info.ExtFunIO.symQs.QsOpp,1) == 0);
    end
    % Muscle activations
    opti.subject_to(a(model_info.ExtFunIO.symQs.MusInvA,end) - a(model_info.ExtFunIO.symQs.MusInvB,1) == 0);
    % Muscle-tendon forces
    opti.subject_to(FTtilde(model_info.ExtFunIO.symQs.MusInvA,end) - FTtilde(model_info.ExtFunIO.symQs.MusInvB,1) == 0);
    % Torque actuator activations
    if ~isempty(model_info.ExtFunIO.symQs.ActInvA)
        opti.subject_to(a_a(model_info.ExtFunIO.symQs.ActInvA,end) - a_a(model_info.ExtFunIO.symQs.ActInvB,1) == 0);
    end
    if ~isempty(model_info.ExtFunIO.symQs.ActOpp)
        opti.subject_to(a_a(model_info.ExtFunIO.symQs.ActOpp,end) + a_a(model_info.ExtFunIO.symQs.ActOpp,1) == 0);
    end
    % Symmetry constraints for synergies
    if (S.subject.synergies)
        opti.subject_to(SynH_r(:,end) - SynH_l(:,1) == 0);
        opti.subject_to(SynH_l(:,end) - SynH_r(:,1) == 0);
    end
else
    opti.subject_to(Qs(model_info.ExtFunIO.symQs.QsFullGC,end) - Qs(model_info.ExtFunIO.symQs.QsFullGC,1) == 0);
    opti.subject_to(Qdots(:,end) - Qdots(:,1) == 0);
    % Muscle activations
    opti.subject_to(a(:,end) - a(:,1) == 0);
    % Muscle-tendon forces
    opti.subject_to(FTtilde(:,end) - FTtilde(:,1) == 0);
    % Torque actuator activations
    if nq.torqAct > 0
        opti.subject_to(a_a(:,end) - a_a(:,1) == 0);
    end
    if (S.subject.synergies)
        opti.subject_to(SynH_r(:,end) - SynH_r(:,1) == 0);
        opti.subject_to(SynH_l(:,end) - SynH_l(:,1) == 0);
    end

end
% Average speed
% Provide expression for the distance traveled
Qs_nsc = Qs.*(scaling.Qs'*ones(1,N+1));
dist_trav_tot = Qs_nsc(model_info.ExtFunIO.jointi.base_forward,end) - ...
    Qs_nsc(model_info.ExtFunIO.jointi.base_forward,1);
vel_aver_tot = dist_trav_tot/tf;
opti.subject_to(vel_aver_tot - S.misc.forward_velocity == 0)

% optional constraints
if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    % lower bound on distance traveled
    if ~isempty(S.bounds.dist_trav.lower)
        opti.subject_to(dist_trav_tot >= S.bounds.dist_trav.lower/2)
    end
    % upper bound on step length for right foot (left foot is already symmetric
    if ~isempty(S.bounds.SLL.upper) || ~isempty(S.bounds.SLR.upper)
        [step_length_r,~] = f_casadi.f_getStepLength(Qs_nsc(:,1),Qs_nsc(:,end));
        if ~isempty(S.bounds.SLL.upper)
            opti.subject_to(step_length_r <= S.bounds.SLL.upper)
        elseif ~isempty(S.bounds.SLR.upper)
            opti.subject_to(step_length_r <= S.bounds.SLR.upper)
        end
    end

elseif strcmp(S.misc.gaitmotion_type,'FullGaitCycle')
    % lower bound on distance traveled
    if ~isempty(S.bounds.dist_trav.lower)
        opti.subject_to(dist_trav_tot >= S.bounds.dist_trav.lower)
    end
    % upper bound on step length for left and/or right foot
    if ~isempty(S.bounds.SLL.upper) || ~isempty(S.bounds.SLR.upper)
        [step_length_r,step_length_l] = f_casadi.f_getStepLength(Qs_nsc(:,1),Qs_nsc(:,end));
        if ~isempty(S.bounds.SLL.upper)
            opti.subject_to(step_length_l <= S.bounds.SLL.upper)
        end
        if ~isempty(S.bounds.SLR.upper)
            opti.subject_to(step_length_r <= S.bounds.SLR.upper)
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scale cost function
Jall_sc = sum(Jall)/dist_trav_tot;
opti.minimize(Jall_sc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp(['...OCP formulation done. Time elapsed ' num2str(toc(t0),'%.2f') ' s'])
disp(' ')

%%

if ~S.post_process.load_prev_opti_vars
    % Create NLP solver
    options = S.solver.nlpsol_options;
    options.ipopt = S.solver.ipopt_options;
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy           = 'adaptive';
    options.ipopt.max_iter              = S.solver.max_iter;
    options.ipopt.linear_solver         = S.solver.linear_solver;
    options.ipopt.tol                   = 1*10^(-S.solver.tol_ipopt);
    options.ipopt.constr_viol_tol       = 1*10^(-S.solver.constr_viol_tol_ipopt);
    opti.solver('ipopt', options);
    % timer
    
    disp('Starting NLP solver...')
    disp(' ')
    t0s = tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve problem
    % Opti does not use bounds on variables but constraints. This function
    % adjusts for that.
    [w_opt,stats] = solve_NLPSOL(opti,options);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(' ')
    disp(['...Exit NLP solver. Time elapsed ' num2str(toc(t0s),'%.2f') ' s'])
    disp(' ')
    disp(' ')
    % Create setup
    setup.tolerance.ipopt = S.solver.tol_ipopt;
    setup.bounds = bounds;
    setup.bounds_nsc = bounds_nsc;
    setup.scaling = scaling;
    setup.guess = guess;
    
    Outname = fullfile(S.misc.save_folder,[S.misc.result_filename '.mat']);
    save(Outname,'w_opt','stats','setup','model_info','S');

else % S.post_process.load_prev_opti_vars = true
    
    % Advanced feature, for debugging only: load w_opt and reconstruct R before rerunning the post-processing.
    Outname = fullfile(S.misc.save_folder,[S.misc.result_filename '.mat']);
    disp(['Loading vector with optimization variables from previous solution: ' Outname])
    clear 'S'
    load(Outname,'w_opt','stats','setup','model_info','R','S');
    scaling = setup.scaling;
    if exist('R','var')
        S = R.S;
    end
    clear 'R'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Essential post processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read from the vector with optimization results
disp('Retrieving solution...')
disp(' ')

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
if nq.torqAct > 0
    a_a_opt = reshape(w_opt(starti:starti+nq.torqAct*(N+1)-1),nq.torqAct,N+1)';
    starti = starti + nq.torqAct*(N+1);
    a_a_col_opt = reshape(w_opt(starti:starti+nq.torqAct*(d*N)-1),nq.torqAct,d*N)';
    starti = starti + nq.torqAct*(d*N);
end
vA_opt = reshape(w_opt(starti:starti+NMuscle*N-1),NMuscle,N)';
starti = starti + NMuscle*N;
if nq.torqAct > 0
    e_a_opt = reshape(w_opt(starti:starti+nq.torqAct*N-1),nq.torqAct,N)';
    starti = starti + nq.torqAct*N;
end
dFTtilde_col_opt=reshape(w_opt(starti:starti+NMuscle*(d*N)-1),NMuscle,d*N)';
starti = starti + NMuscle*(d*N);
qdotdot_col_opt =reshape(w_opt(starti:starti+nq.all*(d*N)-1),nq.all,(d*N))';
starti = starti + nq.all*(d*N);
if (S.subject.synergies)       
    if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
        SynH_r_opt = reshape(w_opt(starti:starti+S.subject.NSyn_r*(N+1)-1),S.subject.NSyn_r,N+1)';
        starti = starti + S.subject.NSyn_r*(N+1);
        SynH_l_opt = reshape(w_opt(starti:starti+S.subject.NSyn_l*(N+1)-1),S.subject.NSyn_l,N+1)';
        starti = starti + S.subject.NSyn_l*(N+1);
        SynW_r_opt = reshape(w_opt(starti:starti+NMuscle/2*S.subject.NSyn_r-1),NMuscle/2,S.subject.NSyn_r)';
        starti = starti + NMuscle/2*S.subject.NSyn_r;
        SynW_l_opt = SynW_r_opt;
    elseif strcmp(S.misc.gaitmotion_type,'FullGaitCycle')
        SynH_r_opt = reshape(w_opt(starti:starti+S.subject.NSyn_r*(N+1)-1),S.subject.NSyn_r,N+1)';
        starti = starti + S.subject.NSyn_r*(N+1);
        SynH_l_opt = reshape(w_opt(starti:starti+S.subject.NSyn_l*(N+1)-1),S.subject.NSyn_l,N+1)';
        starti = starti + S.subject.NSyn_l*(N+1);
        SynW_r_opt = reshape(w_opt(starti:starti+NMuscle/2*S.subject.NSyn_r-1),NMuscle/2,S.subject.NSyn_r)';
        starti = starti + NMuscle/2*S.subject.NSyn_r;
        SynW_l_opt = reshape(w_opt(starti:starti+NMuscle/2*S.subject.NSyn_l-1),NMuscle/2,S.subject.NSyn_l)';
        starti = starti + NMuscle/2*S.subject.NSyn_l;
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
if nq.torqAct > 0
    a_a_mesh_col_opt=zeros(N*(d+1)+1,nq.torqAct);
    a_a_mesh_col_opt(1:(d+1):end,:)= a_a_opt;
end
for k=1:N
    rangei = k*(d+1)-(d-1):k*(d+1);
    rangebi = (k-1)*d+1:k*d;
    a_mesh_col_opt(rangei,:) = a_col_opt(rangebi,:);
    FTtilde_mesh_col_opt(rangei,:) = FTtilde_col_opt(rangebi,:);
    Qs_mesh_col_opt(rangei,:) = Qs_col_opt(rangebi,:);
    Qdots_mesh_col_opt(rangei,:) = Qdots_col_opt(rangebi,:);
    if nq.torqAct > 0
        a_a_mesh_col_opt(rangei,:) = a_a_col_opt(rangebi,:);
    end
end

%% Unscale results
% States at mesh points
% Qs (1:N-1)
q_opt_unsc.rad = Qs_opt(1:end-1,:).*repmat(...
    scaling.Qs,size(Qs_opt(1:end-1,:),1),1);
% Convert in degrees
q_opt_unsc.deg = q_opt_unsc.rad;
q_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations) ...
    = q_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations).*180/pi;
% Qs (1:N)
q_opt_unsc_all.rad = Qs_opt.*repmat(scaling.Qs,size(Qs_opt,1),1);
% Convert in degrees
q_opt_unsc_all.deg = q_opt_unsc_all.rad;
q_opt_unsc_all.deg(:,model_info.ExtFunIO.jointi.rotations) ...
    = q_opt_unsc_all.deg(:,model_info.ExtFunIO.jointi.rotations).*180/pi;
% Qdots (1:N-1)
qdot_opt_unsc.rad = Qdots_opt(1:end-1,:).*repmat(...
    scaling.Qdots,size(Qdots_opt(1:end-1,:),1),1);
% Convert in degrees
qdot_opt_unsc.deg = qdot_opt_unsc.rad;
qdot_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations) ...
    = qdot_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations).*180/pi;
% Qdots (1:N)
qdot_opt_unsc_all.rad =Qdots_opt.*repmat(scaling.Qdots,size(Qdots_opt,1),1);
% Muscle activations (1:N-1)
a_opt_unsc = a_opt(1:end-1,:).*repmat(...
    scaling.a,size(a_opt(1:end-1,:),1),size(a_opt,2));
% Muscle-tendon forces (1:N-1)
FTtilde_opt_unsc = FTtilde_opt(1:end-1,:).*repmat(...
    scaling.FTtilde,size(FTtilde_opt(1:end-1,:),1),1);
if nq.torqAct > 0
    % Torque actuator activations (1:N-1)
    a_a_opt_unsc = a_a_opt(1:end-1,:);
    % Torque actuator activations (1:N)
    a_a_opt_unsc_all = a_a_opt;
end
if (S.subject.synergies)
    % Muscle synergies (1:N-1) at mesh points
    SynH_r_opt_unsc = SynH_r_opt(1:end-1,:).*repmat(...
        scaling.a,size(SynH_r_opt(1:end-1,:),1),size(SynH_r_opt,2)); % same scaling as a
    SynH_l_opt_unsc = SynH_l_opt(1:end-1,:).*repmat(...
        scaling.a,size(SynH_l_opt(1:end-1,:),1),size(SynH_l_opt,2)); % same scaling as a
end

% Controls at mesh points
% Time derivative of muscle activations (states)
vA_opt_unsc = vA_opt.*repmat(scaling.vA,size(vA_opt,1),size(vA_opt,2));

% Torque actuator excitations
if nq.torqAct > 0
    e_a_opt_unsc = e_a_opt;
end
% States at collocation points
% Qs
q_col_opt_unsc.rad = Qs_col_opt.*repmat(scaling.Qs,size(Qs_col_opt,1),1);
% Convert in degrees
q_col_opt_unsc.deg = q_col_opt_unsc.rad;
q_col_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations) ...
    = q_col_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations).*180/pi;
% Qdots
qdot_col_opt_unsc.rad = Qdots_col_opt.*repmat(...
    scaling.Qdots,size(Qdots_col_opt,1),1);
% Convert in degrees
qdot_col_opt_unsc.deg = qdot_col_opt_unsc.rad;
qdot_col_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations) ...
    = qdot_col_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations).*180/pi;
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
qdotdot_col_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations) ...
    = qdotdot_col_opt_unsc.deg(:,model_info.ExtFunIO.jointi.rotations).*180/pi;
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

if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    t_mesh_GC = [tgrid(1:end-1),tgrid(1:end)+tgrid(end)];
else
    t_mesh_GC = tgrid;
end

%% Assert average speed
dist_trav_opt = q_opt_unsc_all.rad(end,model_info.ExtFunIO.jointi.base_forward) - ...
    q_opt_unsc_all.rad(1,model_info.ExtFunIO.jointi.base_forward); % distance traveled
time_elaps_opt = tf_opt; % time elapsed
vel_aver_opt = dist_trav_opt/time_elaps_opt;
% assert_v_tg should be 0
assert_v_tg = abs(vel_aver_opt-S.misc.forward_velocity);
if assert_v_tg > 1*10^(-S.solver.tol_ipopt)
    disp('Issue when reconstructing average speed')
end

%% Decompose optimal cost
J_opt           = 0;
E_cost          = 0;
A_cost          = 0;
Actu_cost       = 0;
Qdotdot_cost    = 0;
Pass_cost       = 0;
vA_cost         = 0;
dFTtilde_cost   = 0;
QdotdotArm_cost = 0;
Syn_cost        = 0;
TrackSyn_cost   = 0;
count           = 1;
h_opt           = tf_opt/N;
for k=1:N
    for j=1:d
        % Get muscle-tendon lengths, velocities, moment arms
        [lMTkj_opt_all,vMTkj_opt_all,~] = ...
            f_casadi.lMT_vMT_dM(q_col_opt_unsc.rad(count,:),qdot_col_opt_unsc.rad(count,:));
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
            MuscleMass',pctsts,full(Fiso_opt_all),model_info.mass,S.metabolicE.tanh_b);
        else
            error('No energy model selected');
        end
        e_tot_opt_all = full(e_tot_all)';
        
        % passive moments
        Tau_passkj = full(f_casadi.AllPassiveTorques_cost(q_col_opt_unsc.rad(count,:),qdot_col_opt_unsc.rad(count,:)));

        % objective function
        J_opt = J_opt + 1/(dist_trav_opt)*(...
            W.E*B(j+1)          *(f_casadi.J_muscles_exp(e_tot_opt_all,W.E_exp))/model_info.mass*h_opt + ...
            W.a*B(j+1)          *(f_casadi.J_muscles_exp(a_col_opt(count,:), W.a_exp))*h_opt + ...
            W.q_dotdot*B(j+1)   *(f_casadi.J_not_arms_dof(qdotdot_col_opt(count,model_info.ExtFunIO.jointi.noarmsi)))*h_opt + ...
            W.pass_torq*B(j+1)  *(f_casadi.J_lim_torq(Tau_passkj))*h_opt + ... 
            W.slack_ctrl*B(j+1) *(f_casadi.J_muscles(vA_opt(k,:)))*h_opt + ...
            W.slack_ctrl*B(j+1) *(f_casadi.J_muscles(dFTtilde_col_opt(count,:)))*h_opt);
            
        if nq.torqAct > 0
            J_opt = J_opt + 1/(dist_trav_opt)*(W.e_torqAct*B(j+1)      *(f_casadi.J_torq_act(e_a_opt(k,:)))*h_opt);

            Actu_cost = Actu_cost + W.e_torqAct*B(j+1)*(f_casadi.J_torq_act(e_a_opt(k,:)))*h_opt;
        end
        if nq.arms > 0
            J_opt = J_opt + 1/(dist_trav_opt)*(W.slack_ctrl*B(j+1) *(f_casadi.J_arms_dof(qdotdot_col_opt(count,model_info.ExtFunIO.jointi.armsi)))*h_opt);

            QdotdotArm_cost = QdotdotArm_cost + W.slack_ctrl*B(j+1)*...
                (f_casadi.J_arms_dof(qdotdot_col_opt(count,model_info.ExtFunIO.jointi.armsi)))*h_opt;
        end

        if (S.subject.synergies)
            syn_constr_k_r = a_opt(k,idx_m_r) - SynH_r_opt(k,:)*SynW_r_opt;
            syn_constr_k_l = a_opt(k,idx_m_l) - SynH_l_opt(k,:)*SynW_l_opt;
            J_opt = J_opt + 1/(dist_trav_opt)*W.SynConstr * B(j+1) *(f_casadi.J_muscles([syn_constr_k_r,syn_constr_k_l]))*h_opt;

            Syn_cost = Syn_cost + W.SynConstr * B(j+1) *(f_casadi.J_muscles([syn_constr_k_r,syn_constr_k_l]))*h_opt;
        end

        E_cost = E_cost + W.E*B(j+1)*...
            (f_casadi.J_muscles_exp(e_tot_opt_all,W.E_exp))/model_info.mass*h_opt;
        A_cost = A_cost + W.a*B(j+1)*...
            (f_casadi.J_muscles(a_col_opt(count,:)))*h_opt;      
        Qdotdot_cost = Qdotdot_cost + W.q_dotdot*B(j+1)*...
            (f_casadi.J_not_arms_dof(qdotdot_col_opt(count,model_info.ExtFunIO.jointi.noarmsi)))*h_opt;
        Pass_cost = Pass_cost + W.pass_torq*B(j+1)*...
            (f_casadi.J_lim_torq(Tau_passkj))*h_opt;
        vA_cost = vA_cost + W.slack_ctrl*B(j+1)*...
            (f_casadi.J_muscles(vA_opt(k,:)))*h_opt;
        dFTtilde_cost = dFTtilde_cost + W.slack_ctrl*B(j+1)*...
            (f_casadi.J_muscles(dFTtilde_col_opt(count,:)))*h_opt;
        count = count + 1;
    end
end

if (S.subject.synergies) && (S.subject.TrackSynW)
    TrackSyn_cost = W.TrackSynW*f_casadi.TrackSynW(SynW_r_opt', SynW_l_opt');
    J_opt = J_opt + N/dist_trav_opt*TrackSyn_cost;
end

J_optf = full(J_opt);
E_costf = full(E_cost);
A_costf = full(A_cost);
Arm_costf = full(Actu_cost);
Qdotdot_costf = full(Qdotdot_cost);
Pass_costf = full(Pass_cost);
vA_costf = full(vA_cost);
dFTtilde_costf = full(dFTtilde_cost);
QdotdotArm_costf = full(QdotdotArm_cost);
Syn_costf = full(Syn_cost);
TrackSyn_costf = full(TrackSyn_cost);

contributionCost.absoluteValues = 1/(dist_trav_opt)*[E_costf,A_costf,...
    Arm_costf,Qdotdot_costf,Pass_costf,vA_costf,dFTtilde_costf,...
    QdotdotArm_costf,Syn_costf,N*TrackSyn_costf];
contributionCost.relativeValues = 1/(dist_trav_opt)*[E_costf,A_costf,...
    Arm_costf,Qdotdot_costf,Pass_costf,vA_costf,dFTtilde_costf,...
    QdotdotArm_costf,Syn_costf,N*TrackSyn_costf]./J_optf*100;
contributionCost.relativeValuesRound2 = ...
    round(contributionCost.relativeValues,2);
contributionCost.labels = {'metabolic energy','muscle activation',...
    'actuator excitation','joint accelerations','limit torques','dadt','dFdt',...
    'arm accelerations','synergy constraints','synergy weights tracking'};

% assertCost should be 0
assertCost = abs(J_optf - 1/(dist_trav_opt)*(E_costf+A_costf + Arm_costf + ...
    Qdotdot_costf + Pass_costf + vA_costf + dFTtilde_costf + QdotdotArm_costf + Syn_costf + N*TrackSyn_costf ));
assertCost2 = abs(stats.iterations.obj(end) - J_optf);

if assertCost > 1*10^(-S.solver.tol_ipopt)
    disp('Issue when reconstructing optimal cost wrt sum of terms')
    disp(['   Difference = ' num2str(assertCost)])
end
if assertCost2 > 1*10^(-S.solver.tol_ipopt)
    disp('Issue when reconstructing optimal cost wrt stats')
    disp(['   Difference = ' num2str(assertCost2)])
end

% % Test collocation function
% [coll_eq_constr_opt,coll_ineq_constr1_opt,coll_ineq_constr2_opt,coll_ineq_constr3_opt,...
%     coll_ineq_constr4_opt,coll_ineq_constr5_opt,coll_ineq_constr6_opt,Jall_opt] = f_coll_map(tf_opt,...
%     a_opt(1:end-1,:)', a_col_opt', FTtilde_opt(1:end-1,:)', FTtilde_col_opt', Qs_opt(1:end-1,:)', ...
%     Qs_col_opt', Qdots_opt(1:end-1,:)', Qdots_col_opt', a_a_opt(1:end-1,:)', a_a_col_opt', ...
%     vA_opt', e_a_opt', dFTtilde_col_opt', qdotdot_col_opt');
% 
% Jall_sc_opt = full(sum(Jall_opt)/dist_trav_opt);
% assertCost3 = abs(stats.iterations.obj(end) - Jall_sc_opt);
% 
% if assertCost3 > 1*10^(-S.solver.tol_ipopt)
%     disp('Issue when reconstructing optimal cost wrt stats')
%     disp(['   Difference = ' num2str(assertCost3)])
% end

%% Reconstruct full gait cycle

% joint accelerations controls on mesh points (2:N)
qddot_opt_unsc.deg = qdotdot_col_opt_unsc.deg(d:d:end,:);
qddot_opt_unsc.rad = qdotdot_col_opt_unsc.rad(d:d:end,:);

if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    % Use symmetry to reconstruct 2nd half of the gait cycle

    idx_2nd_half_GC = N+1:2*N;

    % Qs
    q_opt_unsc.deg = [q_opt_unsc.deg(1:end,:); q_opt_unsc.deg(1:end,:)];
    q_opt_unsc.deg(idx_2nd_half_GC,model_info.ExtFunIO.symQs.QsInvA) = q_opt_unsc.deg(idx_2nd_half_GC,model_info.ExtFunIO.symQs.QsInvB);
    q_opt_unsc.deg(idx_2nd_half_GC,model_info.ExtFunIO.symQs.QsOpp) = -q_opt_unsc.deg(idx_2nd_half_GC,model_info.ExtFunIO.symQs.QsOpp);
    q_opt_unsc.deg(idx_2nd_half_GC,model_info.ExtFunIO.jointi.base_forward) = q_opt_unsc.deg(idx_2nd_half_GC,model_info.ExtFunIO.jointi.base_forward) + dist_trav_opt;

    q_opt_unsc.rad = q_opt_unsc.deg;
    q_opt_unsc.rad(:,model_info.ExtFunIO.jointi.rotations) = q_opt_unsc.rad(:,model_info.ExtFunIO.jointi.rotations).*pi/180;

    dist_trav_opt = dist_trav_opt*2;

    % Qdots
    qdot_opt_unsc.deg = [qdot_opt_unsc.deg(1:end,:); qdot_opt_unsc.deg(1:end,:)];
    qdot_opt_unsc.deg(idx_2nd_half_GC,model_info.ExtFunIO.symQs.QdotsInvA) = qdot_opt_unsc.deg(idx_2nd_half_GC,model_info.ExtFunIO.symQs.QdotsInvB);
    qdot_opt_unsc.deg(idx_2nd_half_GC,model_info.ExtFunIO.symQs.QsOpp) = -qdot_opt_unsc.deg(idx_2nd_half_GC,model_info.ExtFunIO.symQs.QsOpp);

    qdot_opt_unsc.rad = qdot_opt_unsc.deg;
    qdot_opt_unsc.rad(:,model_info.ExtFunIO.jointi.rotations) = qdot_opt_unsc.rad(:,model_info.ExtFunIO.jointi.rotations).*pi/180;

    % Qdotdots
    qddot_opt_unsc.deg = [qddot_opt_unsc.deg; qddot_opt_unsc.deg(1:end,:)];
    qddot_opt_unsc.deg(idx_2nd_half_GC,model_info.ExtFunIO.symQs.QdotsInvA) = qddot_opt_unsc.deg(idx_2nd_half_GC,model_info.ExtFunIO.symQs.QdotsInvB);
    qddot_opt_unsc.deg(idx_2nd_half_GC,model_info.ExtFunIO.symQs.QsOpp) = -qddot_opt_unsc.deg(idx_2nd_half_GC,model_info.ExtFunIO.symQs.QsOpp);

    qddot_opt_unsc.rad = qddot_opt_unsc.deg;
    qddot_opt_unsc.rad(:,model_info.ExtFunIO.jointi.rotations) = qddot_opt_unsc.rad(:,model_info.ExtFunIO.jointi.rotations).*pi/180;

    % Muscle activations
    a_opt_unsc = [a_opt_unsc(1:end,:); a_opt_unsc(1:end,:)];
    a_opt_unsc(idx_2nd_half_GC,model_info.ExtFunIO.symQs.MusInvA) = a_opt_unsc(idx_2nd_half_GC,model_info.ExtFunIO.symQs.MusInvB);

    % Time derivatives of muscle activations
    vA_opt_unsc = [vA_opt_unsc(1:end,:); vA_opt_unsc(1:end,:)];
    vA_opt_unsc(idx_2nd_half_GC,model_info.ExtFunIO.symQs.MusInvA) = vA_opt_unsc(idx_2nd_half_GC,model_info.ExtFunIO.symQs.MusInvB);

    % Muscle-tendon forces
    FTtilde_opt_unsc = [FTtilde_opt_unsc(1:end,:); FTtilde_opt_unsc(1:end,:)];
    FTtilde_opt_unsc(idx_2nd_half_GC,model_info.ExtFunIO.symQs.MusInvA) = FTtilde_opt_unsc(idx_2nd_half_GC,model_info.ExtFunIO.symQs.MusInvB);

    % Time derivative of muscle-tendon force
    dFTtilde_opt_unsc = [dFTtilde_opt_unsc; dFTtilde_opt_unsc(1:end,:)];
    dFTtilde_opt_unsc(idx_2nd_half_GC,model_info.ExtFunIO.symQs.MusInvA) = dFTtilde_opt_unsc(idx_2nd_half_GC,model_info.ExtFunIO.symQs.MusInvB);

    if nq.torqAct > 0
        % Torque actuator activations
        a_a_opt_unsc = [a_a_opt_unsc(1:end,:); a_a_opt_unsc(1:end,:)];
        a_a_opt_unsc(idx_2nd_half_GC,model_info.ExtFunIO.symQs.ActInvA) = a_a_opt_unsc(idx_2nd_half_GC,model_info.ExtFunIO.symQs.ActInvB);
        a_a_opt_unsc(idx_2nd_half_GC,model_info.ExtFunIO.symQs.ActOpp) = -a_a_opt_unsc(idx_2nd_half_GC,model_info.ExtFunIO.symQs.ActOpp);

        % Torque actuator excitations
        e_a_opt_unsc = [e_a_opt_unsc(1:end,:); e_a_opt_unsc(1:end,:)];
        e_a_opt_unsc(idx_2nd_half_GC,model_info.ExtFunIO.symQs.ActInvA) = e_a_opt_unsc(idx_2nd_half_GC,model_info.ExtFunIO.symQs.ActInvB);
        e_a_opt_unsc(idx_2nd_half_GC,model_info.ExtFunIO.symQs.ActOpp) = -e_a_opt_unsc(idx_2nd_half_GC,model_info.ExtFunIO.symQs.ActOpp);

    end
    
    % Synergy activations
    if (S.subject.synergies)
        SynH_r_opt_unsc_half = SynH_r_opt_unsc;
        SynH_l_opt_unsc_half = SynH_l_opt_unsc;
        SynH_r_opt_unsc = [SynH_r_opt_unsc_half; SynH_l_opt_unsc_half]; % mesh points
        SynH_l_opt_unsc = [SynH_l_opt_unsc_half; SynH_r_opt_unsc_half];
    end

end

% express slack controls on mesh points 1:N to be consistent
qddot_opt_unsc.deg = [qddot_opt_unsc.deg(end,:); qddot_opt_unsc.deg(1:end-1,:)];
qddot_opt_unsc.rad = [qddot_opt_unsc.rad(end,:); qddot_opt_unsc.rad(1:end-1,:)];
dFTtilde_opt_unsc = [dFTtilde_opt_unsc(end,:); dFTtilde_opt_unsc(1:end-1,:)];

%% Gait cycle starts at right side initial contact

% Ground reaction forces at mesh points (1:N-1)
Foutk_opt                   = zeros(size(q_opt_unsc.rad,1),F.nnz_out);
for i = 1:size(q_opt_unsc.rad,1)
    % Create zero input vector for external function
    F_ext_input = zeros(model_info.ExtFunIO.input.nInputs,1);
    % Assign Qs
    F_ext_input(model_info.ExtFunIO.input.Qs.all,1) = q_opt_unsc.rad(i,:);
    % Assign Qdots
    F_ext_input(model_info.ExtFunIO.input.Qdots.all,1) = qdot_opt_unsc.rad(i,:);
    % Assign Qdotdots (A)
    F_ext_input(model_info.ExtFunIO.input.Qdotdots.all,1) = qddot_opt_unsc.rad(i,:);

    % Evaluate external function
    res = F(F_ext_input);
    Foutk_opt(i,:) = full(res);
end
GRFk_opt = Foutk_opt(:,[model_info.ExtFunIO.GRFs.right_total model_info.ExtFunIO.GRFs.left_total]);


[idx_GC,idx_GC_base_forward_offset,HS1,HS_threshold] = getStancePhaseSimulation(GRFk_opt,model_info.mass/3);

Qs_GC = q_opt_unsc.deg(idx_GC,:);
Qdots_GC = qdot_opt_unsc.deg(idx_GC,:);
Qdotdots_GC = qddot_opt_unsc.deg(idx_GC,:);
Acts_GC = a_opt_unsc(idx_GC,:);
dActs_GC = vA_opt_unsc(idx_GC,:);
FTtilde_GC = FTtilde_opt_unsc(idx_GC,:);
dFTtilde_GC = dFTtilde_opt_unsc(idx_GC,:);
if nq.torqAct > 0
    a_a_GC = a_a_opt_unsc(idx_GC,:);
    e_a_GC = e_a_opt_unsc(idx_GC,:);
end

if(S.subject.synergies)
    SynH_r_GC = SynH_r_opt_unsc(idx_GC,:);
    SynH_l_GC = SynH_l_opt_unsc(idx_GC,:);
end

% adjust forward position to be continuous and start at 0
Qs_GC(idx_GC_base_forward_offset,model_info.ExtFunIO.jointi.base_forward) = Qs_GC(idx_GC_base_forward_offset,model_info.ExtFunIO.jointi.base_forward) + dist_trav_opt;
Qs_GC(:,model_info.ExtFunIO.jointi.base_forward) = Qs_GC(:,model_info.ExtFunIO.jointi.base_forward) - Qs_GC(1,model_info.ExtFunIO.jointi.base_forward);

%% Unscale actuator torques
if nq.torqAct > 0
    T_a_GC = zeros(size(a_a_GC));
    for i=1:nq.torqAct
        T_a_GC(:,i) = a_a_GC(:,i).*scaling.ActuatorTorque(i);
    end
end

%% Save the results
% Structure Results_all
R.S = S;
R.objective = contributionCost;
R.time.mesh = tgrid;
R.time.coll = tgrid_ext;
R.time.mesh_GC = t_mesh_GC;
R.colheaders.coordinates = model_info.ExtFunIO.coord_names.all;
R.colheaders.muscles = model_info.muscle_info.muscle_names;
R.colheaders.objective = contributionCost.labels;
R.kinematics.Qs = Qs_GC;
R.kinematics.Qdots = Qdots_GC;
R.kinematics.Qddots = Qdotdots_GC;
R.muscles.a = Acts_GC;
R.muscles.da = dActs_GC;
if (S.subject.synergies)
    R.muscles.SynH_r = SynH_r_GC;   
    R.muscles.SynH_l = SynH_l_GC;    
    if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
        R.muscles.SynW_r = SynW_r_opt;
        R.muscles.SynW_l = SynW_r_opt;
    elseif strcmp(S.misc.gaitmotion_type,'FullGaitCycle')    
        R.muscles.SynW_r = SynW_r_opt;
        R.muscles.SynW_l = SynW_l_opt;
    end
end
R.muscles.FTtilde = FTtilde_GC;
R.muscles.dFTtilde = dFTtilde_GC;
if nq.torqAct > 0
    R.torque_actuators.a = a_a_GC;
    R.torque_actuators.e = e_a_GC;
    R.torque_actuators.T = T_a_GC;
else
    R.torque_actuators.a = [];
    R.torque_actuators.e = [];
    R.torque_actuators.T = [];
end
R.ground_reaction.threshold = HS_threshold;
R.ground_reaction.initial_contact_side = HS1;
R.ground_reaction.idx_GC = idx_GC;
R.spatiotemp.dist_trav = dist_trav_opt;

% save results
Outname = fullfile(S.misc.save_folder,[S.misc.result_filename '.mat']);
disp(['Saving results as: ' Outname])
save(Outname,'w_opt','stats','setup','R','model_info');


end

