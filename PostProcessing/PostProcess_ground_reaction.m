function [R] = PostProcess_ground_reaction(S,model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcess_ground_reaction
%   This function calculates the ground reaction forces and moments by 
%   evaluating the external function for the optimal kinematics and adds 
%   the results to the struct with results.
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
% Original date: 13/May/2022
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

% Ground Reaction Forces
R.GRFs = Foutk_opt(:,[model_info.ExtFunIO.GRFs.right_foot,model_info.ExtFunIO.GRFs.left_foot]);
R.colheaders.GRF = {'fore_aft_r','vertical_r','lateral_r','fore_aft_l','vertical_l','lateral_l'};

% Ground Reaction Moments
R.GRMs = Foutk_opt(:,[model_info.ExtFunIO.GRMs.right_total,model_info.ExtFunIO.GRMs.left_total]);

% Center of Pressure
COP = zeros(size(R.GRFs));
for i = 1:N
    if R.GRFs(i,2) > 20
        COP(i,1) = R.GRMs(i,3)./R.GRFs(i,2);
        COP(i,3) = -R.GRMs(i,1)./R.GRFs(i,2);
    end
    if R.GRFs(i,5) > 20
        COP(i,4) = R.GRMs(i,6)./R.GRFs(i,5);
        COP(i,6) = -R.GRMs(i,4)./R.GRFs(i,5);
    end
end
R.COP = COP;



