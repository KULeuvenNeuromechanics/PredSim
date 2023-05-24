function [R] = PostProcess_ground_reaction(model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcess_ground_reaction
%   This function calculates the ground reaction forces and moments by 
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
% Original date: 13/May/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

N = size(R.kinematics.Qs,1);

import casadi.*

F  = external('F',replace(fullfile(R.S.misc.subject_path,R.S.misc.external_function),'\','/'));

Foutk_opt = zeros(N,F.nnz_out);

for i = 1:N
    % Create zero input vector for external function
    F_ext_input = zeros(model_info.ExtFunIO.input.nInputs,1);
    % Assign Qs
    F_ext_input(model_info.ExtFunIO.input.Qs.all,1) = R.kinematics.Qs_rad(i,:);
    % Assign Qdots
    F_ext_input(model_info.ExtFunIO.input.Qdots.all,1) = R.kinematics.Qdots_rad(i,:);
    % Assign Qdotdots (A)
    F_ext_input(model_info.ExtFunIO.input.Qdotdots.all,1) = R.kinematics.Qdotdots_rad(i,:);

    % Evaluate external function
    res = F(F_ext_input);

    Foutk_opt(i,:) = full(res);
end

% Ground Reaction Forces
R.ground_reaction.GRF_r = Foutk_opt(:,model_info.ExtFunIO.GRFs.right_total);
R.ground_reaction.GRF_l = Foutk_opt(:,model_info.ExtFunIO.GRFs.left_total);
R.colheaders.GRF = {'fore_aft','vertical','lateral'};


% Ground Reaction Moments
R.ground_reaction.GRM_r = Foutk_opt(:,model_info.ExtFunIO.GRMs.right_total);
R.ground_reaction.GRM_l = Foutk_opt(:,model_info.ExtFunIO.GRMs.left_total);

% Center of Pressure
COP_r = zeros(size(R.ground_reaction.GRF_r));
COP_l = COP_r;

idx_stance_r = find(R.ground_reaction.GRF_r(:,2) > R.ground_reaction.threshold);
idx_stance_l = find(R.ground_reaction.GRF_l(:,2) > R.ground_reaction.threshold);

COP_r(idx_stance_r,1) = R.ground_reaction.GRM_r(idx_stance_r,3)...
    ./R.ground_reaction.GRF_r(idx_stance_r,2);
COP_r(idx_stance_r,3) = -R.ground_reaction.GRM_r(idx_stance_r,1)...
    ./R.ground_reaction.GRF_r(idx_stance_r,2);

COP_l(idx_stance_l,1) = R.ground_reaction.GRM_l(idx_stance_l,3)...
    ./R.ground_reaction.GRF_l(idx_stance_l,2);
COP_l(idx_stance_l,3) = -R.ground_reaction.GRM_l(idx_stance_l,1)...
    ./R.ground_reaction.GRF_l(idx_stance_l,2);

R.ground_reaction.COP_r = COP_r;
R.ground_reaction.COP_l = COP_l;

R.ground_reaction.idx_stance_r = idx_stance_r;
R.ground_reaction.idx_stance_l = idx_stance_l;



