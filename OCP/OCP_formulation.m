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
% Original date: Januari-May/2022
%
% Last edit by: 
% Last edit date: 
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
F  = external('F',S.misc.external_function);
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

if strcmp(S.subject.IG_selection,'quasi-random')
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
J           = 0; % Initialize cost function
eq_constr   = {}; % Initialize equality constraint vector
ineq_constr1   = {}; % Initialize inequality constraint vector
ineq_constr2   = {}; % Initialize inequality constraint vector
ineq_constr3   = {}; % Initialize inequality constraint vector
ineq_constr4   = {}; % Initialize inequality constraint vector
ineq_constr5   = {}; % Initialize inequality constraint vector
ineq_constr6   = {}; % Initialize inequality constraint vector
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

    % Add contribution to the cost function
    J = J + ...
        W.E          * B(j+1) *(f_casadi.J_muscles_exp(e_totj,W.E_exp))/model_info.mass*h + ...
        W.a          * B(j+1) *(f_casadi.J_muscles(akj(:,j+1)'))*h + ...
        W.q_dotdot   * B(j+1) *(f_casadi.J_not_arms_dof(Aj(model_info.ExtFunIO.jointi.noarmsi,j)))*h + ...
        W.pass_torq  * B(j+1) *(f_casadi.J_lim_torq(Tau_passj_cost))*h + ...
        W.slack_ctrl * B(j+1) *(f_casadi.J_muscles(vAk))*h + ...
        W.slack_ctrl * B(j+1) *(f_casadi.J_muscles(dFTtildej(:,j)))*h;

    if nq.torqAct > 0
        J = J + W.e_arm      * B(j+1) *(f_casadi.J_torq_act(e_ak))*h;
    end
    if nq.arms > 0
        J = J + W.slack_ctrl * B(j+1) *(f_casadi.J_arms_dof(Aj(model_info.ExtFunIO.jointi.armsi,j)))*h;
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Call external function (run inverse dynamics)
    [Tj] = F([QsQdotskj_nsc(:,j+1);Aj_nsc(:,j)]);

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
        
        % total coordinate torque equals inverse dynamics torque
        eq_constr{end+1} = (Tj(i,1) - Ti)./scaling.Moments(i);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Activation dynamics (implicit formulation)
    act1 = vAk_nsc + akj(:,j+1)./(ones(size(akj(:,j+1),1),1)*tdeact);
    act2 = vAk_nsc + akj(:,j+1)./(ones(size(akj(:,j+1),1),1)*tact);
    ineq_constr1{end+1} = act1;
    ineq_constr2{end+1} = act2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contraction dynamics (implicit formulation)
    eq_constr{end+1} = Hilldiffj;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constraints to prevent parts of the skeleton to penetrate each other.
    % Origins calcaneus (transv plane) at minimum 9 cm from each other.
    if ~isempty(model_info.ExtFunIO.origin.calcn_r) &&  ~isempty(model_info.ExtFunIO.origin.calcn_l)
        Qconstr = f_casadi.J_nn_2(Tj(model_info.ExtFunIO.origin.calcn_r([1 3]),1) - ...
            Tj(model_info.ExtFunIO.origin.calcn_l([1 3]),1));
        ineq_constr3{end+1} = Qconstr;
    end
    % Constraint to prevent the arms to penetrate the skeleton
    % Origins femurs and ipsilateral hands (transv plane) at minimum
    % 18 cm from each other.
    if ~isempty(model_info.ExtFunIO.origin.femur_r) &&  ~isempty(model_info.ExtFunIO.origin.hand_r)
        Qconstr = f_casadi.J_nn_2(Tj(model_info.ExtFunIO.origin.femur_r([1 3]),1) - ...
            Tj(model_info.ExtFunIO.origin.hand_r([1 3]),1));
        ineq_constr4{end+1} = Qconstr;
    end
    if ~isempty(model_info.ExtFunIO.origin.femur_l) &&  ~isempty(model_info.ExtFunIO.origin.hand_l)
        Qconstr = f_casadi.J_nn_2(Tj(model_info.ExtFunIO.origin.femur_l([1 3]),1) - ...
            Tj(model_info.ExtFunIO.origin.hand_l([1 3]),1));
        ineq_constr4{end+1} = Qconstr;
    end
    % Origins tibia (transv plane) at minimum 11 cm from each other.
    if ~isempty(model_info.ExtFunIO.origin.tibia_r) &&  ~isempty(model_info.ExtFunIO.origin.tibia_l)
        Qconstr = f_casadi.J_nn_2(Tj(model_info.ExtFunIO.origin.tibia_r([1 3]),1) - ...
            Tj(model_info.ExtFunIO.origin.tibia_l([1 3]),1));
        ineq_constr5{end+1} = Qconstr;
    end
    % Origins toes (transv plane) at minimum 10 cm from each other.
    if ~isempty(model_info.ExtFunIO.origin.toes_r) &&  ~isempty(model_info.ExtFunIO.origin.toes_l)
        Qconstr = f_casadi.J_nn_2(Tj(model_info.ExtFunIO.origin.toes_r([1 3]),1) - ...
            Tj(model_info.ExtFunIO.origin.toes_l([1 3]),1));
        ineq_constr6{end+1} = Qconstr;
    end

end % End loop over collocation points

eq_constr = vertcat(eq_constr{:});
ineq_constr1 = vertcat(ineq_constr1{:});
ineq_constr2 = vertcat(ineq_constr2{:});
ineq_constr3 = vertcat(ineq_constr3{:});
ineq_constr4 = vertcat(ineq_constr4{:});
ineq_constr5 = vertcat(ineq_constr5{:});
ineq_constr6 = vertcat(ineq_constr6{:});


% Casadi function to get constraints and objective
coll_input_vars_def = {tfk,ak,aj,FTtildek,FTtildej,Qsk,Qsj,Qdotsk,Qdotsj,vAk,dFTtildej,Aj};
if nq.torqAct > 0
    coll_input_vars_def = [coll_input_vars_def,{a_ak,a_aj,e_ak}];
end
f_coll = Function('f_coll',coll_input_vars_def,...
    {eq_constr,ineq_constr1,ineq_constr2,ineq_constr3,...
    ineq_constr4,ineq_constr5,ineq_constr6,J});

% assign function to multiple cores
f_coll_map = f_coll.map(N,S.solver.parallel_mode,S.solver.N_threads);

% evaluate function with opti variables
coll_input_vars_eval = {tf,a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col,...
    Qs(:,1:end-1), Qs_col, Qdots(:,1:end-1), Qdots_col, vA, dFTtilde_col, A_col};
if nq.torqAct > 0
    coll_input_vars_eval = [coll_input_vars_eval, {a_a(:,1:end-1), a_a_col, e_a}];
end
[coll_eq_constr,coll_ineq_constr1,coll_ineq_constr2,coll_ineq_constr3,...
    coll_ineq_constr4,coll_ineq_constr5,coll_ineq_constr6,Jall] = f_coll_map(coll_input_vars_eval{:});

% equality constraints
opti.subject_to(coll_eq_constr == 0);

% inequality constraints (logical indexing not possible in MX arrays)
opti.subject_to(coll_ineq_constr1(:) >= 0);
opti.subject_to(coll_ineq_constr2(:) <= 1/tact);
if ~isempty(coll_ineq_constr3)
    opti.subject_to(S.bounds.calcn_dist.lower.^2 < coll_ineq_constr3(:) < 4);
else
    disp('   Minimal distance between calcanei not constrained. To do so, please use "calcn_r" and "calcn_l" as body names in the OpenSim model.')
end
if ~isempty(coll_ineq_constr4)
    opti.subject_to(S.bounds.femur_hand_dist.lower.^2 < coll_ineq_constr4(:) < 4);
end
if isempty(model_info.ExtFunIO.origin.femur_r) || isempty(model_info.ExtFunIO.origin.hand_r)
    disp('   Minimal distance between right arm and body not constrained. To do so, please use "femur_r" and "hand_r" as body names in the OpenSim model.')
end
if isempty(model_info.ExtFunIO.origin.femur_l) || isempty(model_info.ExtFunIO.origin.hand_l)
    disp('   Minimal distance between left arm and body not constrained. To do so, please use "femur_l" and "hand_l" as body names in the OpenSim model.')
end
if ~isempty(coll_ineq_constr5)
    opti.subject_to(S.bounds.tibia_dist.lower.^2 < coll_ineq_constr5(:) < 4);
else
    disp('   Minimal distance between tibias not constrained. To do so, please use "tibia_r" and "tibia_l" as body names in the OpenSim model.')
end
if ~isempty(coll_ineq_constr6)
    opti.subject_to(S.bounds.toes_dist.lower.^2 < coll_ineq_constr6(:) < 4);
else
    disp('   Minimal distance between toes not constrained. To do so, please use "toes_r" and "toes_l" as body names in the OpenSim model.')
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
    opti.subject_to(Qs(model_info.ExtFunIO.symQs.QsOpp,end) + Qs(model_info.ExtFunIO.symQs.QsOpp,1) == 0);
    opti.subject_to(Qdots(model_info.ExtFunIO.symQs.QsOpp,end) + Qdots(model_info.ExtFunIO.symQs.QsOpp,1) == 0);
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

end
% Average speed
% Provide expression for the distance traveled
Qs_nsc = Qs.*(scaling.Qs'*ones(1,N+1));
dist_trav_tot = Qs_nsc(model_info.ExtFunIO.jointi.base_forward,end) - ...
    Qs_nsc(model_info.ExtFunIO.jointi.base_forward,1);
vel_aver_tot = dist_trav_tot/tf;
opti.subject_to(vel_aver_tot - S.subject.v_pelvis_x_trgt == 0)

% optional constraints
if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    % lower bound on distance traveled
    if ~isempty(S.bounds.dist_trav.lower)
        opti.subject_to(dist_trav_tot >= S.bounds.dist_trav.lower/2)
    end
    % upper bound on step length for right foot (left foot is already symmetric
    if ~isempty(S.bounds.SLL.upper) || ~isempty(S.bounds.SLR.upper)
        if ~isempty(model_info.ExtFunIO.origin.calcn_r) &&  ~isempty(model_info.ExtFunIO.origin.calcn_l)
            [step_length_r,~] = f_casadi.f_getStepLength(Qs_nsc(:,1),Qs_nsc(:,end));
            if ~isempty(S.bounds.SLL.upper)
                opti.subject_to(step_length_r <= S.bounds.SLL.upper)
            elseif ~isempty(S.bounds.SLR.upper)
                opti.subject_to(step_length_r <= S.bounds.SLR.upper)
            end
        else
            disp('   Unable to constrain step length. To do so, please use "calcn_r" and "calcn_l" as body names in the OpenSim model.')
        end
    end

elseif strcmp(S.misc.gaitmotion_type,'FullGaitCycle')
    % lower bound on distance traveled
    if ~isempty(S.bounds.dist_trav.lower)
        opti.subject_to(dist_trav_tot >= S.bounds.dist_trav.lower)
    end
    % upper bound on step length for left and/or right foot
    if ~isempty(S.bounds.SLL.upper) || ~isempty(S.bounds.SLR.upper)
        if ~isempty(model_info.ExtFunIO.origin.calcn_r) &&  ~isempty(model_info.ExtFunIO.origin.calcn_l)
            [step_length_r,step_length_l] = f_casadi.f_getStepLength(Qs_nsc(:,1),Qs_nsc(:,end));
            if ~isempty(S.bounds.SLL.upper)
                opti.subject_to(step_length_l <= S.bounds.SLL.upper)
            end
            if ~isempty(S.bounds.SLR.upper)
                opti.subject_to(step_length_r <= S.bounds.SLR.upper)
            end
        else
            disp('   Unable to constrain step length. To do so, please use "calcn_r" and "calcn_l" as body names in the OpenSim model.')
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale cost function
Jall_sc = sum(Jall)/dist_trav_tot;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp(['...OCP formulation done. Time elapsed ' num2str(toc(t0),'%.2f') ' s'])
disp(' ')

%%

if ~S.post_process.load_prev_opti_vars
    % Create NLP solver
    opti.minimize(Jall_sc);
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy           = 'adaptive';
    options.ipopt.max_iter              = S.solver.max_iter;
    options.ipopt.linear_solver         = S.solver.linear_solver;
    options.ipopt.tol                   = 1*10^(-S.solver.tol_ipopt);
    options.ipopt.constr_viol_tol       = 1*10^(-S.solver.tol_ipopt);
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
    
    Outname = fullfile(S.subject.save_folder,[S.post_process.result_filename '.mat']);
    save(Outname,'w_opt','stats','setup','model_info','S');

else % S.post_process.load_prev_opti_vars = true
    
    % Advanced feature, for debugging only: load w_opt and reconstruct R before rerunning the post-processing.
    Outname = fullfile(S.subject.save_folder,[S.post_process.result_filename '.mat']);
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
assert_v_tg = abs(vel_aver_opt-S.subject.v_pelvis_x_trgt);
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
            W.a*B(j+1)          *(f_casadi.J_muscles(a_col_opt(count,:)))*h_opt + ...
            W.q_dotdot*B(j+1)   *(f_casadi.J_not_arms_dof(qdotdot_col_opt(count,model_info.ExtFunIO.jointi.noarmsi)))*h_opt + ...
            W.pass_torq*B(j+1)  *(f_casadi.J_lim_torq(Tau_passkj))*h_opt + ... 
            W.slack_ctrl*B(j+1) *(f_casadi.J_muscles(vA_opt(k,:)))*h_opt + ...
            W.slack_ctrl*B(j+1) *(f_casadi.J_muscles(dFTtilde_col_opt(count,:)))*h_opt);
            
        if nq.torqAct > 0
            J_opt = J_opt + 1/(dist_trav_opt)*(W.e_arm*B(j+1)      *(f_casadi.J_torq_act(e_a_opt(k,:)))*h_opt);

            Actu_cost = Actu_cost + W.e_arm*B(j+1)*(f_casadi.J_arms_dof(e_a_opt(k,:)))*h_opt;
        end
        if nq.arms > 0
            J_opt = J_opt + 1/(dist_trav_opt)*(W.slack_ctrl*B(j+1) *(f_casadi.J_arms_dof(qdotdot_col_opt(count,model_info.ExtFunIO.jointi.armsi)))*h_opt);

            QdotdotArm_cost = QdotdotArm_cost + W.slack_ctrl*B(j+1)*...
                (f_casadi.J_arms_dof(qdotdot_col_opt(count,model_info.ExtFunIO.jointi.armsi)))*h_opt;
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
J_optf = full(J_opt);
E_costf = full(E_cost);
A_costf = full(A_cost);
Arm_costf = full(Actu_cost);
Qdotdot_costf = full(Qdotdot_cost);
Pass_costf = full(Pass_cost);
vA_costf = full(vA_cost);
dFTtilde_costf = full(dFTtilde_cost);
QdotdotArm_costf = full(QdotdotArm_cost);

contributionCost.absoluteValues = 1/(dist_trav_opt)*[E_costf,A_costf,...
    Arm_costf,Qdotdot_costf,Pass_costf,vA_costf,dFTtilde_costf,...
    QdotdotArm_costf];
contributionCost.relativeValues = 1/(dist_trav_opt)*[E_costf,A_costf,...
    Arm_costf,Qdotdot_costf,Pass_costf,vA_costf,dFTtilde_costf,...
    QdotdotArm_costf]./J_optf*100;
contributionCost.relativeValuesRound2 = ...
    round(contributionCost.relativeValues,2);
contributionCost.labels = {'metabolic energy','muscle activation',...
    'actuator excitation','joint accelerations','limit torques','dadt','dFdt',...
    'arm accelerations'};

% assertCost should be 0
assertCost = abs(J_optf - 1/(dist_trav_opt)*(E_costf+A_costf + Arm_costf + ...
    Qdotdot_costf + Pass_costf + vA_costf + dFTtilde_costf + QdotdotArm_costf));

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

end

% express slack controls on mesh points 1:N to be consistent
qddot_opt_unsc.deg = [qddot_opt_unsc.deg(end,:); qddot_opt_unsc.deg(1:end-1,:)];
qddot_opt_unsc.rad = [qddot_opt_unsc.rad(end,:); qddot_opt_unsc.rad(1:end-1,:)];
dFTtilde_opt_unsc = [dFTtilde_opt_unsc(end,:); dFTtilde_opt_unsc(1:end-1,:)];

%% Gait cycle starts at right side initial contact

% Ground reaction forces at mesh points (1:N-1)
Xk_Qs_Qdots_opt             = zeros(size(q_opt_unsc.rad,1),2*nq.all);
Xk_Qs_Qdots_opt(:,1:2:end)  = q_opt_unsc.rad(1:end,:);
Xk_Qs_Qdots_opt(:,2:2:end)  = qdot_opt_unsc.rad(1:end,:);
Xk_Qdotdots_opt             = qddot_opt_unsc.rad(1:end,:);
Foutk_opt                   = zeros(size(q_opt_unsc.rad,1),F.nnz_out);
for i = 1:size(q_opt_unsc.rad,1)
    % ID moments
    [res] = F([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)']);
    Foutk_opt(i,:) = full(res);
end
GRFk_opt = Foutk_opt(:,[model_info.ExtFunIO.GRFs.right_foot model_info.ExtFunIO.GRFs.left_foot]);


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


% save results
Outname = fullfile(S.subject.save_folder,[S.post_process.result_filename '.mat']);
disp(['Saving results as: ' Outname])
save(Outname,'w_opt','stats','setup','R','model_info');


end

