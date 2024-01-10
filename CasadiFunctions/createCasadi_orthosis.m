function [f_orthosis_mesh_k, f_orthosis_mesh_all] = createCasadi_orthosis(S,model_info)
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
%
% Original author: Lars D'Hondt
% Original date: 5/January/2024
% --------------------------------------------------------------------------

import casadi.*

n_coord = model_info.ExtFunIO.jointi.nq.all;

% external function
F = external('F',fullfile(S.misc.subject_path, S.misc.external_function));
F_k = F;
F_all = F.map(S.solver.N_meshes,S.solver.parallel_mode,S.solver.N_threads);


for N = [1,S.solver.N_meshes]


%% Define variables
% all inputs
q_SX = SX.sym('q',n_coord,N);
qdot_SX = SX.sym('qdot',n_coord,N);
qddot_SX = SX.sym('qddot',n_coord,N);
act_SX = SX.sym('a',model_info.muscle_info.NMuscle,N);
fromExtFun_SX = SX.sym('fromExtFun',F.size1_out(0),N); % GRFs, point kinematics

% all outputs
Mcoordk = SX(n_coord,N);
toExtFun_SX = SX(F.size1_in(0),N); % bodyforces, bodymoments


%% loop over all selected orthoses
if ~isempty(S.orthosis)
    for i=1:length(S.orthosis.settings)

        orthosis_i = S.orthosis.settings{i}.object;
        Nmesh = orthosis_i.getNmesh();
        if Nmesh==N

            % Get casadi Function of this orthosis
            f_orthosis_i = orthosis_i.getFunction();
            % allocate input
            arg_i = cell(1,f_orthosis_i.n_in);

            % Use the name of each input to relate it to a variable
            for j=1:f_orthosis_i.n_in
                name_j = f_orthosis_i.name_in(j-1);
                if contains(name_j,'coord_')
                    name_j = replace(name_j,'coord_','');
                    if contains(name_j,'_pos')
                        name_j = replace(name_j,'_pos','');
                        arg_i{j} = q_SX(model_info.ExtFunIO.coordi.(name_j),:);
                    elseif contains(name_j,'_vel')
                        name_j = replace(name_j,'_vel','');
                        arg_i{j} = qdot_SX(model_info.ExtFunIO.coordi.(name_j),:);
                    elseif contains(name_j,'_acc')
                        name_j = replace(name_j,'_acc','');
                        arg_i{j} = qddot_SX(model_info.ExtFunIO.coordi.(name_j),:);
                    end

                elseif contains(name_j,'point_')
                    name_j = replace(name_j,'point_','');
                    if contains(name_j,'_pos')
                        name_j = replace(name_j,'_pos','');
                        arg_i{j} = fromExtFun_SX(model_info.ExtFunIO.position.(name_j),:);
                    elseif contains(name_j,'_vel')
                        name_j = replace(name_j,'_vel','');
                        arg_i{j} = fromExtFun_SX(model_info.ExtFunIO.velocity.(name_j),:);
                    end

                elseif contains(name_j,'GRF_')
                    name_j = replace(name_j,'GRF_','');
                    if strcmp(name_j(end-1:end),'_F')
                        name_j = name_j(1:end-2);
                        arg_i{j} = fromExtFun_SX(model_info.ExtFunIO.GRFs.(name_j),:);
                    end

                elseif contains(name_j,'muscle_')
                    name_j = replace(name_j,'muscle_','');
                    if contains(name_j,'_act')
                        name_j = replace(name_j,'_act','');
                        arg_i{j} = act_SX(strcmp(model_info.muscle_info.muscle_names,name_j),:);
                    end

                end
            end

            % allocate output
            res_i = cell(1,f_orthosis_i.n_out);

            % evaluate casadi Function
            if isempty(arg_i)
                [res_0] = f_orthosis_i();
                for j=1:f_orthosis_i.n_out
                    res_i{j} = res_0.(f_orthosis_i.name_out(j-1));
                end
            else
                [res_i{:}] = f_orthosis_i(arg_i{:});
            end

            % Use the name of each output to relate it to a variable
            for j=1:f_orthosis_i.n_out

                name_j = f_orthosis_i.name_out(j-1);
                if contains(name_j,'CoordForce_')
                    name_j = replace(name_j,'CoordForce_','');
                    Mcoordk(model_info.ExtFunIO.coordi.(name_j)) = Mcoordk(model_info.ExtFunIO.coordi.(name_j)) + res_i{j};

                elseif contains(name_j,'BodyForce_')
                    name_j = replace(name_j,'BodyForce_','');
                    toExtFun_SX(model_info.ExtFunIO.input.Forces.(name_j),:) = toExtFun_SX(model_info.ExtFunIO.input.Forces.(name_j),:) + res_i{j};

                elseif contains(name_j,'BodyMoment_')
                    name_j = replace(name_j,'BodyMoment_','');
                    toExtFun_SX(model_info.ExtFunIO.input.Moments.(name_j),:) = toExtFun_SX(model_info.ExtFunIO.input.Moments.(name_j),:) + res_i{j};

                end
            end
        end
    end
end


%% create casadi Function for combination of orthoses
fun = Function(['f_Orthosis_mesh_',num2str(N),'_wo_ext'],...
    {q_SX, qdot_SX, qddot_SX, act_SX, fromExtFun_SX},...
    {Mcoordk, toExtFun_SX},...
    {'qs','qdots','qddots','act','fromExtFun'},...
    {'M_coord','toExtFun'});

%% Create casadi Function for combination of orthoses, and include external function
% inputs
qs_MX = MX.sym('qs',n_coord,N);
qdots_MX = MX.sym('qdots',n_coord,N);
qddots_MX = MX.sym('qddots',n_coord,N);
act_MX = MX.sym('a',model_info.muscle_info.NMuscle,N);


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
[Mcoord_MX, toExtFun_MX] = fun(qs_MX,qdots_MX,qddots_MX,act_MX,fromExtFun_MX);

% Create function to be used in OCP formulation
if N==1
    f_orthosis_mesh_k = Function('f_orthosis_mesh_k',{qs_MX,qdots_MX,qddots_MX,act_MX},...
        {Mcoord_MX, toExtFun_MX},{'qs','qdots','qddots','act'},{'M_coord','M_body'});
else
    f_orthosis_mesh_all = Function('f_orthosis_mesh_all',{qs_MX,qdots_MX,qddots_MX,act_MX},...
        {Mcoord_MX, toExtFun_MX},{'qs','qdots','qddots','act'},{'M_coord','M_body'});
end

end % for N = [1, Nmeshes]

%% Cleanup
% delete Orthosis objects, because we don't need them anymore and don't want to save it
% note: might want to handle this more cleanly
if ~isempty(S.orthosis)
    for i=1:length(S.orthosis.settings)
%         S.orthosis.settings{i} = rmfield(S.orthosis.settings{i},'object');
        S.orthosis.settings{i}.object.delete();
    end
end

end % end of function