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
[F] = load_external_function(R.S);


QsQdots = zeros(N,2*model_info.ExtFunIO.jointi.nq.all);
QsQdots(:,1:2:end) = R.kinematics.Qs_rad;
QsQdots(:,2:2:end) = R.kinematics.Qdots_rad;

Foutk_opt = zeros(N,F.nnz_out);

for i = 1:N
    % ID moments
    [res] = F([QsQdots(i,:)';R.kinematics.Qddots_rad(i,:)']);
    Foutk_opt(i,:) = full(res);
end

% Ground Reaction Forces
R.colheaders.GRF = {'fore_aft','vertical','lateral'};
GRFs = fieldnames(model_info.ExtFunIO.GRFs);

R.ground_reaction.GRF_r = Foutk_opt(:,model_info.ExtFunIO.GRFs.right_foot);
GRFs = GRFs(~strcmp(GRFs(:),'right_foot'));
R.ground_reaction.GRF_l = Foutk_opt(:,model_info.ExtFunIO.GRFs.left_foot);
GRFs = GRFs(~strcmp(GRFs(:),'left_foot'));

for i=1:length(GRFs)
    R.ground_reaction.(GRFs{i}) = Foutk_opt(:,model_info.ExtFunIO.GRFs.(GRFs{i}));
end

% Stance phase
idx = find(R.ground_reaction.GRF_r(:,2) < 20,1,'first');
idx_stance_r = (1:idx-1)';
trh = min(R.ground_reaction.GRF_r(idx_stance_r,2));
idx_stance_l = find(R.ground_reaction.GRF_l(:,2) > trh);

R.ground_reaction.idx_stance_r = idx_stance_r;
R.ground_reaction.idx_stance_l = idx_stance_l;

if isfield(model_info.ExtFunIO,'GRMs')
    % Ground Reaction Moments
    R.ground_reaction.GRM_r = Foutk_opt(:,model_info.ExtFunIO.GRMs.right_total);
    R.ground_reaction.GRM_l = Foutk_opt(:,model_info.ExtFunIO.GRMs.left_total);
    
    % Center of Pressure
    COP_r = zeros(size(R.ground_reaction.GRF_r));
    COP_l = COP_r;
    
    
    
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
end




