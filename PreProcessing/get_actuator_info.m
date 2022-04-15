function [model_info] = get_actuator_info(S,osim_path,model_info)
% --------------------------------------------------------------------------
% get_actuator_info
%   Read max torque, activation time constant, and coordinate of the ideal
%   torque actuators from the osim file. Actuators should be added as
%   "ActivationCoordinateActuator" objects in opensim. 
%   An example of adding such actuators to an osim file can be found in
%   \Tests\add_actuators_to_osim.m 
%   
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% Original author: Lars D'Hondt
% Original date: 18/March/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------



import org.opensim.modeling.*;

% get all actuators from osim model
model = Model(osim_path);
model.initSystem;

actuators = model.getActuators();
Nact = actuators.getSize();

% prepare struct to write inforation
coord_names = {};
coordi = [];
max_torque = [];
time_constant = [];

% loop over all actuators
for i=1:Nact
    % only use proper type of actuator
    act_type = actuators.get(i-1).getConcreteClassName;
    if strcmp(act_type,'ActivationCoordinateActuator')
        act_i = actuators.get(i-1);
        acta_i = ActivationCoordinateActuator.safeDownCast(act_i);
        % read parameters and add to struct
        coord_name_i = char(acta_i.getCoordinate().getName());
        coord_names{end+1} = coord_name_i;
        max_torque(end+1) = double(acta_i.getOptimalForce());
        time_constant(end+1) = double(acta_i.get_activation_time_constant());
        coordi(end+1) = model_info.ExtFunIO.coordi.(coord_name_i);
    end
end

actuator_info.coord_names = coord_names;
actuator_info.coordi = coordi;
actuator_info.max_torque = max_torque;
actuator_info.time_constant = time_constant;
actuator_info.NActuators = length(actuator_info.coordi);

% place struct with actuator info in model_info
model_info.actuator_info = actuator_info;

% Coordinates actuated by ideal actuators
model_info.ExtFunIO.jointi.torqueActuated = actuator_info.coordi;
model_info.ExtFunIO.jointi.nq.TorqAct = Nact;


