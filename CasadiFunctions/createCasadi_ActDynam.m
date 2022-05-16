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
% Last update:
%   reformulate for generic actuator
% Last edit by: Lars D'Hondt
% Last edit date: 15/April/2022
% --------------------------------------------------------------------------

import casadi.*
time_constants = struct_array_to_double_array(model_info.actuator_info.parameters,'time_constant');

%% Activation dynamics
e = SX.sym('e',model_info.ExtFunIO.jointi.nq.torqAct);
a = SX.sym('a',model_info.ExtFunIO.jointi.nq.torqAct);
dadt = (e-a)./time_constants;

f_ActuatorActivationDynamics = Function('f_ActuatorActivationDynamics',{e,a},{dadt},{'e','a'},{'dadt'});

