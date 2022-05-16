function [R] = PostProcess_external_function(S,model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcess_external_function
%   This function evaluates the external function for the optimal
%   kinematics and adds the results to the struct with results.
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

N = size(R.Qs,1);

import casadi.*
[F] = load_external_function(S);


QsQdots = zeros(N,2*model_info.ExtFunIO.jointi.nq.all);
QsQdots(:,1:2:end) = R.Qs_rad;
QsQdots(:,2:2:end) = R.Qdots_rad;

Foutk_opt = zeros(N,F.nnz_out);

for i = 1:N
    % ID moments
    [res] = F([QsQdots(i,:)';R.Qddots_rad(i,:)']);
    Foutk_opt(i,:) = full(res);
end

R.Tid = Foutk_opt(:,1:model_info.ExtFunIO.jointi.nq.all);
R.colheaders.joints = model_info.ExtFunIO.coord_names.all;

R.GRFs = Foutk_opt(:,[model_info.ExtFunIO.GRFs.right_foot,model_info.ExtFunIO.GRFs.left_foot]);
if size(R.GRFs,2)==6
    R.colheaders.GRF = {'fore_aft_r','vertical_r','lateral_r','fore_aft_l','vertical_l','lateral_l'};
elseif size(R.GRFs,2)==4
    R.colheaders.GRF = {'fore_aft_r','vertical_r','fore_aft_l','vertical_l'};
else
    R.colheaders.GRF = {};
end





end