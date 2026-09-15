function [R] = PostProcess_joint_stiffness(model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcessing_joint_stiffness
%   Function to compute the joint stiffness for all joints
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
%   * struct with joint stiffness field
% 
% Original author: Menthy Denayer
% Original date: 05/June/2026
%
% Last edit by: 
% Last edit date: 
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

%% Define Variables
a = R.muscles.a';
lMtilde = R.muscles.lMtilde';
lT = R.muscles.lT';
lTtilde = lT./[model_info.muscle_info.parameters.lTs]';
lMT = R.muscles.lMT';
vM = R.muscles.vM';
FT = R.muscles.FT';
FM = R.muscles.Fce'+R.muscles.Fpass';
Q = R.kinematics.Qs_rad';
rij = permute(R.muscles.dM, [2 3 1]);

%% Compute Muscle Stiffness
[KT, KM] = f_casadi.f_muscle_tendon_stiffness(a,lMtilde,vM,lT);

KT_full = full(KT);
KM_full = full(KM);

R.joint_stiffness.KT = KT_full';
R.joint_stiffness.KM = KM_full';

%% Compute Derivative of Moment Arm
drdtheta = zeros(size(a,1),size(Q,1),size(Q,2));

for j = 1:size(drdtheta,3)                                                  % loop over all timepoints
    [~,~,~,drdthetaj] =  f_casadi.lMT_vMT_dM(Q(:,j)',[]);                   % use lMT_vMT_dM function to compute muscle moment arm partial derivatives
    drdthetaj_full = full(drdthetaj);
    for i = 1:size(a,1)                                                     % loop over all muscles
        H_i = drdthetaj_full(i:size(a,1):end, :);                           % Nq x Nq Hessian for muscle i: H_i(j,k) = d(dM(i,j))/dtheta_k
        dMdr_diag(i,:) = diag(H_i)';                                        % keep only j = k terms (second order derivatives)
    end
    drdtheta(:,:,j) = dMdr_diag;
end

R.joint_stiffness.drdtheta = permute(drdtheta, [3 1 2]);

%% Compute Joint Stiffness
KMJ = zeros(size(Q,2), size(a,1), size(Q,1));                                
KJ = zeros(size(Q))';                                                          
R.joint_stiffness.KJ = KJ;                                                    
R.joint_stiffness.KMJ = KMJ;

for j = 1:size(KJ,1)
    [KJj, KMJj] = f_casadi.f_joint_stiffness(KM_full(:,j),KT_full(:,j),FT(:,j),FM(:,j),lMT(:,j),...
        lTtilde(:,j),lMtilde(:,j),rij(:,:,j),drdtheta(:,:,j));

    R.joint_stiffness.KJ(j,:) = full(KJj)';
    R.joint_stiffness.KMJ(j,:,:) = full(KMJj);
end