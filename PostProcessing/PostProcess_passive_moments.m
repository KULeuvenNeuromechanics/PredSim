function [R] = PostProcess_passive_moments(model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcess_passive_moments
%   This function computes the different contributions to the passive joint
%   moments.
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
% Original date: 19/May/2022
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


N = size(R.kinematics.Qs,1);

R.kinetics.T_spring = zeros(N,model_info.ExtFunIO.jointi.nq.all);
R.kinetics.T_damping = zeros(N,model_info.ExtFunIO.jointi.nq.all);
R.kinetics.T_limit = zeros(N,model_info.ExtFunIO.jointi.nq.all);

for i=1:N
    Tki = f_casadi.PassiveStiffnessMoments(R.kinematics.Qs_rad(i,:)');
    R.kinetics.T_spring(i,:) = full(Tki)';

    Tdi = f_casadi.PassiveDampingMoments(R.kinematics.Qdots_rad(i,:)');
    R.kinetics.T_damping(i,:) = full(Tdi)';

    Tli = f_casadi.LimitTorques(R.kinematics.Qs_rad(i,:)');
    R.kinetics.T_limit(i,:) = full(Tli)';

end




