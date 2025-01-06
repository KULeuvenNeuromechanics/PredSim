function [] = removeActuators(osim_path)
% --------------------------------------------------------------------------
% removeActuators
%   This functions removes the actuators (ActivationCoordinateActuator) 
%   from a given .osim model.
% 
% INPUT:
%   - osim_path -
%   * path to the OpenSim model file (.osim)
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
        force_j.delete();
    end
end

%% save model
model.finalizeConnections();
model.initSystem();
model.print(osim_path);



end % end of function
