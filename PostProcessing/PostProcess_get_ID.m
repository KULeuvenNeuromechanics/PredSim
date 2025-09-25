function [R] = PostProcess_get_ID(model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcess_get_ID
%   This function calculates the inverse dynamic joint torques/forces by 
%   evaluating the external function for the optimal kinematics and adds 
%   the results to the struct with results.
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
% Original date: 10/May/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

N = size(R.kinematics.Qs,1);

import casadi.*

if R.S.OpenSimADOptions.useSerialisedFunction
    F = Function.load(replace(fullfile(R.S.misc.subject_path,R.S.misc.external_function),'\','/'));
else
    F = external('F',replace(fullfile(R.S.misc.subject_path,R.S.misc.external_function),'\','/'));
end

Foutk_opt = zeros(N,F.nnz_out);

for i = 1:N
    % Create zero input vector for external function
    F_ext_input = zeros(model_info.ExtFunIO.input.nInputs,1);
    % Assign Qs
    F_ext_input(model_info.ExtFunIO.input.Qs.all,1) = R.kinematics.Qs_rad(i,:);
    % Assign Qdots
    F_ext_input(model_info.ExtFunIO.input.Qdots.all,1) = R.kinematics.Qdots_rad(i,:);
    % Assign Qdotdots (A)
    F_ext_input(model_info.ExtFunIO.input.Qdotdots.all,1) = R.kinematics.Qddots_rad(i,:);

    % Add contact forces for kinematic constraints (action-reaction)
    for j=1:model_info.ExtFunIO.jointi.nq.constr
        F_ext_input(model_info.ExtFunIO.input.Forces. ...
            (['osimConstraint_',model_info.osimConstraints{j},'_1']),1) = ...
            R.kinetics.F_kin_constr(i, (1:3)+3*(j-1));
        F_ext_input(model_info.ExtFunIO.input.Forces. ...
            (['osimConstraint_',model_info.osimConstraints{j},'_2']),1) = ...
            - R.kinetics.F_kin_constr(i, (1:3)+3*(j-1));
    end

    % Evaluate external function
    res = F(F_ext_input);

    Foutk_opt(i,:) = full(res);
end

R.kinetics.T_ID = Foutk_opt(:,1:model_info.ExtFunIO.jointi.nq.all);


