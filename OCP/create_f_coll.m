function [f_coll] = create_f_coll(S,model_info,f_casadi,scaling)
% --------------------------------------------------------------------------
% create_f_coll
%   This function creates f_coll, a CasADi Function to evaluate the 
%   constraints and cost function contribution at each collocation point.
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
%   - scaling -
%   * scale factors for all optimisation variables
%
% OUTPUT:
%   - f_coll -
%   * CasADi Function to evaluate the constraints and cost function
%   contribution at each collocation point.
% 
% Original author: Dhruv Gupta and Lars D'Hondt
% Original date: January-May/2022
% --------------------------------------------------------------------------

%% User inputs (typical settings structure)
% settings for optimization
N = S.solver.N_meshes; % number of mesh intervals
W = S.weights; % weights optimization
nq = model_info.ExtFunIO.jointi.nq; % lengths of coordinate subsets
NMuscle = model_info.muscle_info.NMuscle; % Total number of muscles

%% Load external functions
import casadi.*
% The external function performs inverse dynamics through the
% OpenSim/Simbody C++ API. This external function is compiled as a dll from
% which we create a Function instance using CasADi in MATLAB. More details
% about the external function can be found in the documentation.
if S.OpenSimADOptions.useSerialisedFunction
    F = Function.load(fullfile(S.misc.subject_path, S.misc.external_function));
else
    F = external('F',fullfile(S.misc.subject_path, S.misc.external_function));
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

[~,mai] = MomentArmIndices_asym(muscleNames,...
    model_info.muscle_info.muscle_spanning_joint_info);
% calculate total number of joints that each muscle crosses (used later)
sumCross = sum(model_info.muscle_info.muscle_spanning_joint_info);

% Parameters for activation dynamics
tact = model_info.muscle_info.tact; % Activation time constant
tdeact = model_info.muscle_info.tdeact; % Deactivation time constant


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
    SynW_rk         = MX.sym('SynW_rk',length(model_info.muscle_info.idx_right),S.subject.NSyn_r);
    SynW_lk         = MX.sym('SynW_lk',length(model_info.muscle_info.idx_left),S.subject.NSyn_l);
end

% Contact forces for kinematic constraints
if nq.constr > 0
    F_constrj = MX.sym('F_constrj',3*nq.constr,d);
    F_constrj_nsc = F_constrj.*scaling.F_constr;
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
        lMTj,vMTj);
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
            akj(:,j+1),akj(:,j+1),lMtildej,vMj,Fcej,Fpassj,Fisoj);
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
        W.slack_ctrl * B(j+1) *(f_casadi.J_muscles(vAk))*h + ...
        W.slack_ctrl * B(j+1) *(f_casadi.J_muscles(dFTtildej(:,j)))*h;

    if nq.limTorq > 0
        J = J + W.pass_torq  * B(j+1) *(f_casadi.J_lim_torq(Tau_passj_cost))*h;
    end
    if nq.torqAct > 0
        J = J + W.e_torqAct  * B(j+1) *(f_casadi.J_torq_act(e_ak))*h;
    end
    if nq.arms > 0
        J = J + W.slack_ctrl * B(j+1) *(f_casadi.J_arms_dof(Aj(model_info.ExtFunIO.jointi.armsi,j)))*h;
    end
    if nq.constr > 0
        J = J + W.slack_ctrl * B(j+1) *( sumsqr(F_constrj(:,j))/(nq.constr*3) )*h;
    end

    % If muscle synergies: Instead of a - WH = 0 as an equality constraint, have it as a
        % term in the cost function to be minimized (+ inequality constraint)     
    if (S.subject.synergies)        
            syn_constr_k_r = ak(model_info.muscle_info.idx_right) - SynW_rk*SynH_rk;
            syn_constr_k_l = ak(model_info.muscle_info.idx_left) - SynW_lk*SynH_lk;
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
    
    % Add contact forces for kinematic constraints (action-reaction)
    for i=1:nq.constr
        F_ext_input(model_info.ExtFunIO.input.Forces. ...
            (['osimConstraint_',model_info.osimConstraints{i},'_1']),1) = ...
            F_constrj_nsc((1:3)+3*(i-1),j);
        F_ext_input(model_info.ExtFunIO.input.Forces. ...
            (['osimConstraint_',model_info.osimConstraints{i},'_2']),1) = ...
            - F_constrj_nsc((1:3)+3*(i-1),j);
    end

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
%         if ~ismember(i,model_info.ExtFunIO.jointi.floating_base)
            Ti = Ti + Tau_passj(i);
%         end
        
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

    % kinematic constraints
    for i=1:nq.constr
        pos_1 = Tj(model_info.ExtFunIO.position.(['osimConstraint_',model_info.osimConstraints{i},'_1']),1);
        pos_2 = Tj(model_info.ExtFunIO.position.(['osimConstraint_',model_info.osimConstraints{i},'_2']),1);
        eq_constr{end+1} = (pos_2 - pos_1);
    end

end % End loop over collocation points

% Add tracking terms in the cost function if synergy weights are tracked
% Here we select the weights that we want to impose/track 
% (there are no conditions/constraints applied to the other weights)
if (S.subject.synergies) && (S.subject.TrackSynW)
    J_TrackSynW = W.TrackSynW*f_casadi.TrackSynW(SynW_rk, SynW_lk);
    J = J + J_TrackSynW;
else
    J_TrackSynW = 0;
end

% Synergies: a - WH = 0
% Only applied for mesh points
if (S.subject.synergies)
    ineq_constr_syn{end+1} =  [ak(model_info.muscle_info.idx_right);ak(model_info.muscle_info.idx_left)] - [SynW_rk*SynH_rk;SynW_lk*SynH_lk];
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
if nq.constr > 0
    coll_input_vars_def = [coll_input_vars_def,{F_constrj}];
end


f_coll = Function('f_coll',coll_input_vars_def,...
        {eq_constr, ineq_constr_deact, ineq_constr_act,...
        ineq_constr_distance{:}, ineq_constr_syn,J});


if S.OpenSimADOptions.useSerialisedFunction
    try
        f_coll = f_coll.expand();
    catch
    end
end



end
