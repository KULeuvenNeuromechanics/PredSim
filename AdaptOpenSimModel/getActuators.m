function [torq_act] = getActuators(osim_path)
% --------------------------------------------------------------------------
% getActuators
%   This functions reads the parameter values describing the actuators
%   (ActivationCoordinateActuator) of a given .osim model, and returns this
%   information in a structured way.
% 
% INPUT:
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
% OUTPUT:
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
% Original date: 3 January 2025
% --------------------------------------------------------------------------


%% load model
import org.opensim.modeling.*;
model = Model(osim_path);

%% loop over forces
torq_act = [];
i = 1;
for j=1:model.getForceSet().getSize()
    force_j = model.getForceSet().get(j-1);
    if strcmp(force_j.getConcreteClassName(),'ActivationCoordinateActuator')
        actu = ActivationCoordinateActuator.safeDownCast(force_j);

        torq_act(i).coord = char( actu.get_coordinate() );
        torq_act(i).max_torque = actu.getOptimalForce();
        torq_act(i).time_constant = actu.get_activation_time_constant();

        i = i+1;


    end
end




end % end of function
