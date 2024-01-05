function [f_orthosis_mesh_k, f_orthosis_mesh_all] = createCasadi_orthosis(S,model_info)
% --------------------------------------------------------------------------
% createCasadi_orthosis
%   This function creates a casadifunction to calculate the total orthosis
%   torque on every joint.
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
%   - f_orthosis -
%   * casadi function to calculate orthosis torque
% 
% Original author: Lars D'Hondt
% Original date: 14/September/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.*

n_coord = model_info.ExtFunIO.jointi.nq.all;

% external function
F = external('F',fullfile(S.misc.subject_path, S.misc.external_function));


%% Define variables
% all inputs
qk = SX.sym('q',n_coord,1);
qdotk = SX.sym('qdot',n_coord,1);
qddotk = SX.sym('qddot',n_coord,1);
actk = SX.sym('a',model_info.muscle_info.NMuscle,1);
GRFk = {};
pointPosk = {};
pointVelk = {};

% all outputs
Mcoordk = SX(n_coord,1);
Fbodyk = {};
Mbodyk = {};



if ~isempty(S.orthosis)
    for i=1:length(S.orthosis.settings)

        orthosis_i = S.orthosis.settings{i}.object;
        Nmesh = orthosis_i.getNmesh();
        if Nmesh==1

            f_orthosis_i = orthosis_i.getFunction();
            arg_i = cell(1,f_orthosis_i.n_in);

            for j=1:f_orthosis_i.n_in
                name_j = f_orthosis_i.name_in(j-1);
                if contains(name_j,'coord_')
                    name_j = replace(name_j,'coord_','');
                    if contains(name_j,'_pos')
                        name_j = replace(name_j,'_pos','');
                        arg_i{j} = qk(model_info.ExtFunIO.coordi.(name_j));
                    elseif contains(name_j,'_vel')
                        name_j = replace(name_j,'_vel','');
                        arg_i{j} = qdotk(model_info.ExtFunIO.coordi.(name_j));
                    elseif contains(name_j,'_acc')
                        name_j = replace(name_j,'_acc','');
                        arg_i{j} = qddotk(model_info.ExtFunIO.coordi.(name_j));
                    end

                elseif contains(name_j,'point_')
                    name_j = replace(name_j,'point_','');
                    if contains(name_j,'_pos')
                        name_j = replace(name_j,'_pos','');
                        if ~isfield(pointPosk,name_j)
                            pointPosk.(name_j) = SX.sym(name_j,3,1);
                        end
                        arg_i{j} = pointPosk.(name_j);
                    elseif contains(name_j,'_vel')
                        name_j = replace(name_j,'_vel','');
                        if ~isfield(pointVelk,name_j)
                            pointVelk.(name_j) = SX.sym(name_j,3,1);
                        end
                        arg_i{j} = pointPosk.(name_j);
                    end

                elseif contains(name_j,'GRF_')
                    name_j = replace(name_j,'GRF_','');
                    if ~isfield(GRFk,name_j)
                        GRFk.(name_j) = SX.sym(name_j,3,1);
                    end
                    arg_i{j} = GRFk.(name_j);

                elseif contains(name_j,'muscle_')
                    name_j = replace(name_j,'muscle_','');
                    if contains(name_j,'_act')
                        name_j = replace(name_j,'_act','');
                        arg_i{j} = actk(strcmp(model_info.muscle_info.muscle_names,name_j));
                    end

                end
            end

            res_i = cell(1,f_orthosis_i.n_out);

            [res_i{:}] = f_orthosis_i(arg_i{:});


            for j=1:f_orthosis_i.n_out

                name_j = f_orthosis_i.name_in(j-1);
                if contains(name_j,'CoordForce_')
                    name_j = replace(name_j,'CoordForce_','');
                    Mcoordk(model_info.ExtFunIO.coordi.(name_j)) = Mcoordk(model_info.ExtFunIO.coordi.(name_j)) + res_i{j};

                elseif contains(name_j,'BodyForce_')
                    name_j = replace(name_j,'BodyForce_','');
                    if isfield(Fbodyk,name_j)
                        Fbodyk.(name_j) = Fbodyk.(name_j) + res_i{j};
                    else
                        Fbodyk.(name_j) = res_i{j};
                    end

                elseif contains(name_j,'BodyMoment_')
                    name_j = replace(name_j,'BodyMoment_','');
                    if isfield(Mbodyk,name_j)
                        Mbodyk.(name_j) = Mbodyk.(name_j) + res_i{j};
                    else
                        Mbodyk.(name_j) = res_i{j};
                    end
                end
            end
        end
    end
end


GRFk_all = SX(3*length(fieldnames(GRFk)));
GRFk_idx = nan(3*length(fieldnames(GRFk)));


%%

f_orthosis = Function('f_orthosis',{qk,qdotk,GRFk},{Tk},{'qk','qdotk','GRFk'},{'Tk'});



end