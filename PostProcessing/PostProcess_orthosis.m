function [R] = PostProcess_orthosis(model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcessing_orthosis
%   Calculate generalised forces from orthoses.
% 
% INPUT:
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
%   - f_casadi -
%   * Struct containing all casadi functions.
%
%   - R -
%   * struct with simulation results
%
% OUTPUT:
%   - R -
%   * struct with simulation results
% 
% Original author: Lars D'Hondt
% Original date: 22 January 2024
% --------------------------------------------------------------------------


N = R.S.solver.N_meshes;

import casadi.*


qs(R.ground_reaction.idx_GC,:) = R.kinematics.Qs_rad;
qdots(R.ground_reaction.idx_GC,:) = R.kinematics.Qdots_rad;
qddots(R.ground_reaction.idx_GC,:) = R.kinematics.Qddots_rad;
acts(R.ground_reaction.idx_GC,:) = R.muscles.a;

ortArg = {};
if isfield(R,'orthosis')
    if isfield(R.orthosis,'states')
        states(R.ground_reaction.idx_GC,:) = R.orthosis.states;
        states = states(1:N,:);
        ortArg{end+1} = states';
    end
    if isfield(R.orthosis,'controls')
        controls(R.ground_reaction.idx_GC,:) = R.orthosis.controls;
        controls = controls(1:N,:);
        ortArg{end+1} = controls';
    end
end

qs = qs(1:N,:);
qdots = qdots(1:N,:);
qddots = qddots(1:N,:);
acts = acts(1:N,:);

[Mcoord1, toExtFun1] = f_casadi.f_orthosis_mesh_k(qs', qdots', qddots', acts',ortArg{:});

[Mcoord2, toExtFun2] = f_casadi.f_orthosis_mesh_all(qs',qdots',qddots',acts',ortArg{:});


F  = external('F',replace(fullfile(R.S.misc.subject_path,R.S.misc.external_function),'\','/'));


% Create zero input vector for external function
F_ext_input = zeros(model_info.ExtFunIO.input.nInputs,N);
% Assign Qs
F_ext_input(model_info.ExtFunIO.input.Qs.all,:) = qs';
% Assign Qdots
F_ext_input(model_info.ExtFunIO.input.Qdots.all,:) = qdots';
% Assign Qdotdots (A)
F_ext_input(model_info.ExtFunIO.input.Qdotdots.all,:) = qddots';

fromExtFun = full(F(F_ext_input));

F_ext_input = F_ext_input + full(toExtFun1) + full(toExtFun2);

Fout_opt = F(F_ext_input);
Fout_opt = full(Fout_opt(1:model_info.ExtFunIO.jointi.nq.all,:));
Mcoord1 = full(Mcoord1);
Mcoord2 = full(Mcoord2);

tmp.T =  Fout_opt - (Mcoord1+Mcoord2);

tmp2 = construct_full_gait_cycle(model_info, tmp);
T_bio = tmp2.T(:,R.ground_reaction.idx_GC)';

Mcoord = R.kinetics.T_ID - T_bio;

R.orthosis.combined.T_coord = Mcoord;
R.orthosis.combined.T_bio = T_bio;

%%


for i=1:length(f_casadi.separate_orthoses)
    
    f_pp_i = f_casadi.separate_orthoses(i).wrap_pp;

    res_fun = cell(1,f_pp_i.n_out);

    [res_fun{:}] = f_pp_i(qs', qdots', qddots', acts',fromExtFun,ortArg{:});

    for j=1:f_pp_i.n_out
        res_pp.x = full(res_fun{j});
        res_pp = construct_full_gait_cycle(model_info, res_pp);
        res_pp_GC = res_pp.x(:,R.ground_reaction.idx_GC)';

        name_ij = f_pp_i.name_out(j-1);
        R.orthosis.separate{i}.(name_ij) = res_pp_GC;
    end

%     R.orthosis.separate{i} = R_sep_i;
end





end
