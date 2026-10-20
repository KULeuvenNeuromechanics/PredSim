function [AFO] = parametricAFO(init, settings_orthosis)
% --------------------------------------------------------------------------
% parametricAFO
%   Defines an ankle-foot orthosis as a rotational spring on the ankle angle
%   and a rotational spring on the MTP angle. Each spring is defined by a
%   stiffness constant.
% 
%
% INPUT:
%   - init -
%   * struct with information used to initialise the Orthosis object.
% 
%   - settings_orthosis -
%   * struct with information about this orthosis, containing the fields:
%       - function_name = 'parametricAFO'  i.e. name of this function           
%       - ankle_stiffness:  ankle stiffness in Nm/rad
%       - mtp_stiffness:    mtp stiffness in Nm/rad
%       - left_right:       'l' for left or 'r' for right
%   Values are set via S.orthosis.settings{i} in main.m, with i the index
%   of the orthosis.
%
%
% OUTPUT:
%   - AFO -
%   * an object of the class Orthosis
% 
% Original author: Lars D'Hondt
% Original date: 5/January/2024
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

% create Orthosis object
AFO = Orthosis('AFO',init);

% read settings that were passed from main.m
k_ankle = settings_orthosis.ankle_stiffness; % ankle stiffness in Nm/rad
k_mtp = settings_orthosis.mtp_stiffness; % mtp stiffness in Nm/rad
side = settings_orthosis.left_right; % 'l' for left or 'r' for right

% get joint angles
q_ankle = AFO.var_coord(['ankle_angle_',side]); % ankle angle in rad;
q_mtp = AFO.var_coord(['mtp_angle_',side]); % MTP angle in rad;

% calculate moments
T_ankle = -k_ankle*q_ankle;
T_mtp = -k_mtp*q_mtp;

% add calculated moments to Orthosis
AFO.addCoordForce(T_ankle, ['ankle_angle_',side])
AFO.addCoordForce(T_mtp, ['mtp_angle_',side])

end