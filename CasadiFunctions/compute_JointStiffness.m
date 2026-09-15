function [K_J, K_M_J] = compute_JointStiffness(K_M, K_T, FT, FM, lMT, lTtilde, lMtilde, dM, drij_dtheta, model_info, lMo_in,...
    lTs_in, alphao_in, MuscMoAsmp)
% --------------------------------------------------------------------------
% compute_JointStiffness 
%   Function to compute the joint stiffness.
% 
% INPUT:
%   - K_M -
%   * muscle stiffness: partial derivative of total muscle force to muscle
%   length
%
%   - K_T -
%   * tendon stiffness: partial derivative of tendon force to tendon length
%
%   - FT -
%   * tendon force
%
%   - FM -
%   * total muscle force (active + passive)
%
%   - lMT -
%   * muscle-tendon length
%
%   - lTtilde -
%   * normalized tendon length
%
%   - lMtilde -
%   * normalized muscle length
%
%   - dM -
%   * muscle-moment arm
%
%   - drij_dtheta -
%   * muscle-moment arm partial derivative to joint angles
%
%   - lMo_in -
%   * optimal fibre lengths
%
%   - lTs_in -
%   * tendon slack lengths
%
%   - alphao_in -
%   * constant pennation angle
%
%   - MuscMoAsmp -
%   * constant pennation angle
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - K_J -
%   * computed joint stiffness (Nm/rad)
% 
% Original author: Menthy Denayer
% Original date: 08/September/2026
%
% --------------------------------------------------------------------------
% This file is part of PredSim.
% 
% PredSim: A Framework for Rapid Predictive Simulations of Locomotion
% Copyright (c) 2026 KU Leuven
% 
% PredSim is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Affero General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% 
% PredSim is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public 
% License for more details.
% 
% You should have received a copy of the GNU Affero General Public License 
% along with PredSim. If not, see <https://www.gnu.org/licenses/>.
% --------------------------------------------------------------------------

% import casadi
import casadi.*

% define muscle properties
N_muscles = model_info.muscle_info.NMuscle;                                 % number of muscles
lMo = ones(N_muscles,1).*lMo_in';                                           % optimal fibre lengths
lTs = ones(N_muscles,1).*lTs_in';                                           % tendon slack lengths
alphao = ones(N_muscles,1).*alphao_in';                                     % constant pennation angles

%% Define Variables
Njoints = model_info.ExtFunIO.jointi.nq.all;                                % number of joints
lT = lTtilde .* lTs;                                                        % tendon slack lengths (non-normalized)
lM = lMtilde .* lMo;                                                        % muscle fibre lengths (non-normalized)

%% Compute Pennation Angle
if(MuscMoAsmp == 0) % b = cst
    cos_alpha = (lMT-lT)./lM;
else    % alpha = cst = alphao
    cos_alpha = cos(alphao);
end

%% Compute Effective Muscle Stiffness
% Projecting stiffness onto tendon axis
K_Eff = K_M .* cos_alpha.^2 + FM./lM .* (1-cos_alpha.^2);

%% Compute MTU Stiffness
% Springs in series
K_MTU = 1./(1./K_T + 1./K_Eff);

%% Compute Stiffness Muscle-Joint
% K_M_J = K_MTU * rij² + drij/dtheta * F_T

K_MTU_3D = repmat(K_MTU, [1 Njoints]);
FT_3D    = repmat(FT,    [1 Njoints]);

term_elastic = K_MTU_3D .* dM.^2;
term_geometric = drij_dtheta .* FT_3D;

K_M_J = term_elastic + term_geometric;

%% Compute Joint Stiffness
% Sum over all muscle contributions for one joint
% K_J = sum_M K_M_J

K_J = SX.zeros(1,Njoints); % Joint stiffness
for j = 1:Njoints
    K_J(:,j) = sum(K_M_J(:,j), 1);
end

end