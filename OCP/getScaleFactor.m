function [scaling] = getScaleFactor(S,model_info,bounds_nsc)
% --------------------------------------------------------------------------
% getScaleFactor
%   This script provides scaling factors for the optimisation variables.
%   Scale factors are based on the bound with highest absolute value, so
%   optimisation variables remain within the interval [-1,1].
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
%   - bounds_nsc -
%   * boundaries for all optimisation variables (not scaled)
%
% OUTPUT:
%   - scaling -
%   * scale factors for all optimisation variables
% 
% Original author: Lars D'Hondt
% Original date: 6/April/2023
%
% Last edit by:
% Last edit date:
% --------------------------------------------------------------------------


% Qs
scaling.Qs = max(abs(bounds_nsc.Qs.lower),abs(bounds_nsc.Qs.upper));

% Qdots
scaling.Qdots = max(abs(bounds_nsc.Qdots.lower),abs(bounds_nsc.Qdots.upper));

% Qdotdots
scaling.Qdotdots = max(abs(bounds_nsc.Qdotdots.lower),abs(bounds_nsc.Qdotdots.upper));

% Torque actuator 
scaling.ActuatorTorque = struct_array_to_double_array(model_info.actuator_info.parameters,'max_torque');

% Time derivative of muscle activations
scaling.vA = 100;

% Muscle activations
scaling.a = 1;

% Torque actuator activations
scaling.a_a = 1;

% Torque actuator excitations
scaling.e_a = 1;

% Time derivative of muscle-tendon forces
scaling.dFTtilde = 100;

% Muscle-tendon forces
scaling.FTtilde = max(abs(bounds_nsc.FTtilde.lower),abs(bounds_nsc.FTtilde.upper)); 


end
