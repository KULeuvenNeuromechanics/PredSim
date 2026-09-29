function [f_muscle_tendon_stiffness,f_joint_stiffness] = createCasadi_JointStiffness(S,model_info)
% --------------------------------------------------------------------------
% createCasadi_JointStiffness
%   Function to create Casadi functions for joint stiffness.
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - f_muscle_tendon_stiffness -
%   * function to compute the muscle & tendon stiffness
% 
% Original author: Ines Vandekerckhove, Tom Buurke & Dhruv Gupta, KU Leuven
% Original date: 30-11-2021 
%
% Last edit by: Menthy Denayer
% Last edit date: 01/September/2026 : Removed drdtheta function as included in lMT_vMT_dM Casadi function
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


%% Define Variabkes
import casadi.*
N_muscles = model_info.muscle_info.NMuscle;
N_joints = model_info.ExtFunIO.jointi.nq.all;

%% Muscle & Tendon Stiffness
a          = SX.sym('a',N_muscles);                                         % Muscle activations
lMtilde_in = SX.sym('lMtilde_in',N_muscles);                                % Muscle fibre lengths
vM         = SX.sym('vM',N_muscles);                                        % Muscle fibre velocities
lT         = SX.sym('lT',N_muscles);                                        % Tendon length
KM         = SX(N_muscles,1);                                               % Muscle stiffness
KT         = SX(N_muscles,1);                                               % Tendon stiffness

% Parameters of force-length-velocity curves
load('Ftparam.mat','Ftparam');
load('Fvparam.mat','Fvparam');
load('Fpparam.mat','Fpparam');
load('Faparam.mat','Faparam');

% Function to get muscle & tendons stiffness
for m = 1:N_muscles
    [KT(m), KM(m)] = ForceEquilibrium_dFtildeState_all_tendon(a(m),lMtilde_in(m),vM(m),lT(m),...
        model_info.muscle_info.parameters(m).FMo,model_info.muscle_info.parameters(m).lMo,...
        model_info.muscle_info.parameters(m).lTs,model_info.muscle_info.parameters(m).vMmax,...
        Ftparam,Fvparam,Fpparam,Faparam,...
        model_info.muscle_info.parameters(m).muscle_pass_stiff_shift,...
        model_info.muscle_info.parameters(m).muscle_pass_stiff_scale,model_info.muscle_info.parameters(m).muscle_strength,...
        model_info.muscle_info.parameters(m).tendon_stiff);
end
f_muscle_tendon_stiffness = ...
    Function('f_muscle_tendon_stiffness',{a,lMtilde_in,vM,lT},{KT, KM},...
    {'a','lMtilde','vM','lT'},...
    {'KT','KM'});


%% Joint Stiffness
KM_in       = SX.sym('KM_in',N_muscles);                                    % Muscle stiffness
KT_in       = SX.sym('KT_in',N_muscles);                                    % Tendon stiffness
drdtheta_in = SX.sym('drdtheta_in',N_muscles,N_joints);                     % Derivative of moment arm
FM          = SX.sym('FM',N_muscles);                                       % Muscle force
FT          = SX.sym('FT',N_muscles);                                       % Muscle force
lMT         = SX.sym('lMT',N_muscles);                                      % Muscle-tendon length
lTtilde_in  = SX.sym('lTtilde_in',N_muscles);                               % Tendon length normalized
dM          = SX.sym('dM',N_muscles,N_joints);                              % Muscle moment arms

[KJ,KMJ] = compute_JointStiffness(KM_in, KT_in, FT, FM, lMT, lTtilde_in, lMtilde_in, dM, drdtheta_in, model_info, ...
    [model_info.muscle_info.parameters.lMo],[model_info.muscle_info.parameters.lTs],...
    [model_info.muscle_info.parameters.alphao], S.misc.constant_pennation_angle);

f_joint_stiffness = ...
    Function('f_joint_stiffness',{KM_in,KT_in,FT,FM,lMT,lTtilde_in,lMtilde_in,dM,drdtheta_in},{KJ,KMJ},...
    {'KM_in','KT_in','FT','FM','lMT','lTtilde_in','lMtilde_in','dM','drdtheta_in'},...
    {'KJ','KMJ'});

end