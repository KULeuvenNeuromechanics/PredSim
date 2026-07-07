function [R] = PostProcess_muscle_excitations(model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcess_muscle_excitations
%   This function computes muscle excitations from time derivative of muscle 
%   activations using Raasch's model.
%   More details in De Groote et al. (2009): DOI: 10.1080/10255840902788587
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


% Activation time constant
tact = ones(1,model_info.muscle_info.NMuscle)*model_info.muscle_info.tact; 

% Deactivation time constant
tdeact = ones(1,model_info.muscle_info.NMuscle)*model_info.muscle_info.tdeact; 

% Compute excitations
R.muscles.e = computeExcitationRaasch(R.muscles.a,R.muscles.da,tdeact,tact);

