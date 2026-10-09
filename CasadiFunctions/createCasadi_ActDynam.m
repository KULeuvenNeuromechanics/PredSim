function [f_ActuatorActivationDynamics] = createCasadi_ActDynam(S,model_info)
% --------------------------------------------------------------------------
% createCasadi_ActDynam 
%   Function to create CasADi functions for activation dynamics for
%   actuators.
% 
% INPUT:
%   - S -
%   * setting structure S
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - f_ActuatorActivationDynamics -
%   * formulation activation dynamics for actuators
% 
% Original author: Tom Buurke
% Original date: 30/November/2021
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

import casadi.*
time_constants = struct_array_to_double_array(model_info.actuator_info.parameters,'time_constant');

%% Activation dynamics
e = SX.sym('e',model_info.ExtFunIO.jointi.nq.torqAct);
a = SX.sym('a',model_info.ExtFunIO.jointi.nq.torqAct);
dadt = (e-a)./time_constants;

f_ActuatorActivationDynamics = Function('f_ActuatorActivationDynamics',{e,a},{dadt},{'e','a'},{'dadt'});

end
