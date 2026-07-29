function [exo] = ankleExoEmgProportional(init, settings_orthosis)
% --------------------------------------------------------------------------
% ankleExoEmgProportional
%   Ankle exoskeleton that applies a torque proportional to soleus
%   activation.
%
%
% INPUT:
%   - init -
%   * struct with information used to initialise the Orthosis object.
% 
%   - settings_orthosis -
%   * struct with information about this orthosis, containing the fields:
%       - function_name = 'ankleExoEmgProportional'  i.e. name of this function
%       - gain:             T = -gain*act
%       - left_right:       'l' for left or 'r' for right
%   Values are set via S.orthosis.settings{i} in main.m, with i the index
%   of the orthosis.
%
%
% OUTPUT:
%   - exo -
%   * an object of the class Orthosis
% 
% Original author: Lars D'Hondt
% Original date: 8/January/2024
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
exo = Orthosis('exo',init);


gain = settings_orthosis.gain;
side = settings_orthosis.left_right; % 'l' for left or 'r' for right

emg = exo.var_muscle(['soleus_',side]);

T_ankle = -(emg-0.05)*gain; % only use activation above lower bound (0.05)
exo.addCoordForce(T_ankle,['ankle_angle_',side]);

end