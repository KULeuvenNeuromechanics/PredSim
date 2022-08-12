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
[F] = load_external_function(R.S);

QsQdots = zeros(N,2*model_info.ExtFunIO.jointi.nq.all);
qs = R.kinematics.Qs_rad;
% qs(:,3) = qs(:,3)+0.03;
% QsQdots(:,1:2:end) = R.kinematics.Qs_rad;
QsQdots(:,1:2:end) = qs;
QsQdots(:,2:2:end) = R.kinematics.Qdots_rad;

Foutk_opt = zeros(N,F.nnz_out);

for i = 1:N
    % ID moments
    [res] = F([QsQdots(i,:)';R.kinematics.Qddots_rad(i,:)']);
    Foutk_opt(i,:) = full(res);
end

R.kinetics.T_ID = Foutk_opt(:,1:model_info.ExtFunIO.jointi.nq.all);


end