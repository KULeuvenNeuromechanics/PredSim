function [] = scaleActuators(osim_path, torque_scale_factor)
% --------------------------------------------------------------------------
% scaleActuators
%   This functions scales the optimal force parameter of the actuators
%   (ActivationCoordinateActuator) of a given .osim model
%
%   Note that if there are actuators on translational coordinates, you may 
%   want to scale those differently. In that case this function needs to be
%   adapted to check whether the coordinate is a rotation or translation
%   and use thae proper scale factor (from 2nd input).
% 
% INPUT:
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%   - torque_scale_factor
%   * scale factor for optimal force of actuator. 
% 
%
% Original author: Lars D'Hondt
% Original date: 6 January 2025
% --------------------------------------------------------------------------


%% load model
import org.opensim.modeling.*;
model = Model(osim_path);

%% loop over forces
for j=1:model.getForceSet().getSize()
    force_j = model.getForceSet().get(j-1);
    if strcmp(force_j.getConcreteClassName(),'ActivationCoordinateActuator')
        actu = ActivationCoordinateActuator.safeDownCast(force_j);

        max_torque = actu.getOptimalForce() * torque_scale_factor;
        actu.setOptimalForce(max_torque);

    end
end

%% save model
model.finalizeConnections();
model.initSystem();
model.print(osim_path);


end % end of function
