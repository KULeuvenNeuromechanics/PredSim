function [f_orthosis_mesh_k, f_orthosis_mesh_all, separate_orthoses] = createCasadi_orthosis(S,model_info)
% --------------------------------------------------------------------------
% createCasadi_orthosis
%   This function creates a casadifunction to calculate the total orthosis
%   forces and moments.
% 
% INPUT:
%   - S -
%   * setting structure S
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
%
% OUTPUT:
%   - f_orthosis_mesh_k -
%   * casadi function to calculate forces and moments from orthoses that
%   act the same for each mesh interval (i.e. passive, pure feedback control)
% 
%   - f_orthosis_mesh_all -
%   * casadi function to calculate forces and moments from orthoses that do
%   not act the same for each mesh interval (i.e. feedforward component),
%   or have interdependence between mesh intervals (i.e. delay).
%
%   - separate_orthoses -
%   * cell array with casadi functions of individual orthoses. Used for
%   post-processing.
%
%
% Original author: Lars D'Hondt
% Original date: 5/January/2024
% --------------------------------------------------------------------------

import casadi.*

n_coord = model_info.ExtFunIO.jointi.nq.all;

separate_orthoses = {};

% external function
F = external('F',fullfile(S.misc.subject_path, S.misc.external_function));
F_k = F;
F_all = F.map(S.solver.N_meshes,S.solver.parallel_mode,S.solver.N_threads);

model_info.ExtFunIO.nOutputs = F.size1_out(0);

for N = [S.solver.N_meshes,1] %changed to first do N=Nmeshes


%% Define variables
% all inputs
q_SX = SX.sym('q',n_coord,N);
qdot_SX = SX.sym('qdot',n_coord,N);
qddot_SX = SX.sym('qddot',n_coord,N);
act_SX = SX.sym('a',model_info.muscle_info.NMuscle,N);
fromExtFun_SX = SX.sym('fromExtFun',F.size1_out(0),N); % GRFs, point kinematics

orthStates_SX = SX.sym('x',length(S.orthosis.stateNames_all),N);
orthControls_SX = SX.sym('u',length(S.orthosis.controlNames_all),N);

% all outputs
Mcoordk_SX = SX(n_coord,N);
toExtFun_SX = SX(F.size1_in(0),N); % bodyforces, bodymoments
stateDyn_SX = SX(S.orthosis.Nstates_all,N);


%% loop over all selected orthoses
if ~isempty(S.orthosis)

    for i=1:length(S.orthosis.settings)

        orthosis_i = S.orthosis.settings{i}.object;
        Nmesh = orthosis_i.getNmesh();
        if Nmesh==N

            % Get casadi Function of this orthosis
            [f_orthosis_i, f_orthosis_pp_i] =...
                orthosis_i.wrapCasadiFunction(model_info.ExtFunIO,model_info.muscle_info.muscle_names,S.orthosis.stateNames_all,S.orthosis.controlNames_all);

            % Add to struct for post-processing
            separate_orthoses(i).wrap = f_orthosis_i;
            separate_orthoses(i).wrap_pp = f_orthosis_pp_i;

            % evaluate function 
            [Mcoordk_i, toExtFun_i, stateDyn_i] = f_orthosis_i(q_SX,qdot_SX,qddot_SX,act_SX,fromExtFun_SX,orthStates_SX,orthControls_SX);
            
            % accumulate outputs
            Mcoordk_SX = Mcoordk_SX + Mcoordk_i;
            toExtFun_SX = toExtFun_SX + toExtFun_i;
            stateDyn_SX = stateDyn_SX + stateDyn_i; % Each orthosis uniquely defines a subset of x_dot, zero otherwise

            % create Casadi function for orthosis state dynamics
            separate_orthoses(i).dynamics = Function(['f_Orthosis_dynamics_',num2str(N)],...
                {q_SX, qdot_SX, qddot_SX, act_SX, fromExtFun_SX,orthStates_SX,orthControls_SX},...
                {stateDyn_i},...
                {'qs','qdots','qddots','act','fromExtFun','orthStates','orthControls'},...
                {'orthStateDyn'});                       
        end
    end
end


%% create casadi Function for combination of orthoses
if Nmesh==N
    fun = Function(['f_Orthosis_mesh_',num2str(N),'_wo_ext'],...
        {q_SX, qdot_SX, qddot_SX, act_SX, fromExtFun_SX,orthStates_SX,orthControls_SX},...
        {Mcoordk_SX, toExtFun_SX,stateDyn_SX},...
        {'qs','qdots','qddots','act','fromExtFun','orthStates','orthControls'},...
        {'M_coord','toExtFun','stateDyn'});
else
    % call function with Nmesh rows in the arguments
    [res_Mcoordk_SX, res_toExtFun_SX, res_stateDyn_SX] = fun(repmat(q_SX,1,Nmesh),repmat(qdot_SX,1,Nmesh),...
        repmat(qddot_SX,1,Nmesh),repmat(act_SX,1,Nmesh),...
        repmat(fromExtFun_SX,1,Nmesh),repmat(orthStates_SX,1,Nmesh),...
        repmat(orthControls_SX,1,Nmesh));
    res_Mcoordk_SX_1 = res_Mcoordk_SX(:,1);
    res_toExtFun_SX_1 = res_toExtFun_SX(:,1);
    res_stateDyn_SX_1 = res_stateDyn_SX(:,1);
    % redefine mapping fun
    fun = Function(['f_Orthosis_mesh_',num2str(N),'_wo_ext'],...
        {q_SX, qdot_SX, qddot_SX, act_SX, fromExtFun_SX,orthStates_SX,orthControls_SX},...
        {res_Mcoordk_SX_1, res_toExtFun_SX_1,res_stateDyn_SX_1},...
        {'qs','qdots','qddots','act','fromExtFun','orthStates','orthControls'},...
        {'M_coord','toExtFun','stateDyn'});

end

%% Create casadi Function for combination of orthoses, and include external function
% inputs
qs_MX = MX.sym('qs',n_coord,N);
qdots_MX = MX.sym('qdots',n_coord,N);
qddots_MX = MX.sym('qddots',n_coord,N);
act_MX = MX.sym('a',model_info.muscle_info.NMuscle,N);

orthStates_MX = MX.sym('x',S.orthosis.Nstates_all,N);
orthControls_MX = MX.sym('u',S.orthosis.Ncontrols_all,N);
stateDyn_MX = MX.sym('x_dot',S.orthosis.Nstates_all,N);

% Create zero (sparse) input vector for external function
F_ext_input = MX(model_info.ExtFunIO.input.nInputs,N);
% Assign Qs
F_ext_input(model_info.ExtFunIO.input.Qs.all,:) = qs_MX;
% Assign Qdots
F_ext_input(model_info.ExtFunIO.input.Qdots.all,:) = qdots_MX;
% Assign Qdotdots (A)
F_ext_input(model_info.ExtFunIO.input.Qdotdots.all,:) = qddots_MX;
% Evaluate external function
if N==1
    fromExtFun_MX = F_k(F_ext_input);
else
    fromExtFun_MX = F_all(F_ext_input);
end

% Evaluate orthosis function
[Mcoord_MX, toExtFun_MX, stateDyn_MX] = fun(qs_MX,qdots_MX,qddots_MX,act_MX,fromExtFun_MX,orthStates_MX,orthControls_MX);

% Create function to be used in OCP formulation
if N==1
    f_orthosis_mesh_k = Function('f_orthosis_mesh_k',{qs_MX,qdots_MX,qddots_MX,act_MX,orthStates_MX,orthControls_MX},...
        {Mcoord_MX, toExtFun_MX, stateDyn_MX},{'qs','qdots','qddots','act','orthStates','orthControls'},{'M_coord','M_body','stateDyn'});
else
    f_orthosis_mesh_all = Function('f_orthosis_mesh_all',{qs_MX,qdots_MX,qddots_MX,act_MX,orthStates_MX,orthControls_MX},...
        {Mcoord_MX, toExtFun_MX, stateDyn_MX},{'qs','qdots','qddots','act','orthStates','orthControls'},{'M_coord','M_body','stateDyn'});
end

end % for N = [1, Nmeshes]

%% Cleanup
%TODO: this section was commented out, removal of object must be postponed
% until after orthosis NLP vars and constraints have been added dynamically
% in OCP_formulation.m

% delete Orthosis objects, because we don't need them anymore and don't want to save it
% note: might want to handle this more cleanly
% if ~isempty(S.orthosis)
%     for i=1:length(S.orthosis.settings)
% %         S.orthosis.settings{i} = rmfield(S.orthosis.settings{i},'object');
%         S.orthosis.settings{i}.object.delete();
%     end
% end


end % end of function