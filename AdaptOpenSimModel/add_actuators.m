function [] = add_actuators(osim_path,torq_act)
% --------------------------------------------------------------------------
% add_actuators
%   This functions adds ActivationCoordinateActuator to selected coordinates 
%   of the .osim model.
% 
% INPUT:
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%   - torq_act -
%   * cell array of structs describing an ideal torque actuator. 
%     Example cell:
%       % coordinate to which te actuator applies torque
%       torq_act(1).coord = ('arm_flex_r');
%       % maximum torque
%       torq_act(1).max_torque = 150;
%       % time constant of the activation dynamics
%       torq_act(1).time_constant = 0.035;
% 
%
% Original author: Lars D'Hondt
% Original date: 27/May/2022
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


%% load model
import org.opensim.modeling.*;
model = Model(osim_path);

%% loop over actuators
for i=1:length(torq_act)
    % get coordinate
    coord = model.getCoordinateSet().get(torq_act(i).coord);
    % construct actuator
    actu = ActivationCoordinateActuator();
    % assign coordinate
    actu.setCoordinate(coord);
    % set name
    actu.setName(['actuator_' torq_act(i).coord]);
    % set optimal torque
    actu.setOptimalForce(torq_act(i).max_torque);
    % set activation time constant
    actu.set_activation_time_constant(torq_act(i).time_constant);
    % default activation = 0
    actu.set_default_activation(0);
    % symmetric torque range
    actu.set_min_control(-1);
    actu.set_max_control(1);
    % add to model
    model.addForce(actu);

    % clear variable names
    clear 'actu'

end

%% save model
model.finalizeConnections();
model.initSystem();
model.print(osim_path);







