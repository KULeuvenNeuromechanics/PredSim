% --------------------------------------------------------------------------
% add_actuators_to_osim_file
%   This script adds ideal torque actuators to selected coordinates on an
%   OpenSim model. 
% 
% 
% Original author: Lars D'Hondt
% Original date: 09/May/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------
clear
clc

import org.opensim.modeling.*;
[pathHere,~,~] = fileparts(mfilename('fullpath'));

% osim_filename = 'Hamner_modified';
osim_filename = 'CP3_T0_scaled_MRI_v7_scaledMT_test';

model = Model([pathHere '\' osim_filename '.osim']);
coords = model.getCoordinateSet();

% shoulder left
coord = coords.get('arm_flex_l');
actu = ActivationCoordinateActuator();
actu.setCoordinate(coord);
actu.setName('actuator_arm_flex_l');
actu.setOptimalForce(150);
actu.set_activation_time_constant(0.035);
actu.set_default_activation(0);
actu.set_min_control(-1);
actu.set_max_control(1);
model.addForce(actu);

coord = coords.get('arm_add_l');
actu = ActivationCoordinateActuator();
actu.setCoordinate(coord);
actu.setName('actuator_arm_add_l');
actu.setOptimalForce(150);
actu.set_activation_time_constant(0.035);
actu.set_default_activation(0);
actu.set_min_control(-1);
actu.set_max_control(1);
model.addForce(actu);

coord = coords.get('arm_rot_l');
actu = ActivationCoordinateActuator();
actu.setCoordinate(coord);
actu.setName('actuator_arm_rot_l');
actu.setOptimalForce(150);
actu.set_activation_time_constant(0.035);
actu.set_default_activation(0);
actu.set_min_control(-1);
actu.set_max_control(1);
model.addForce(actu);

% shoulder right
coord = coords.get('arm_flex_r');
actu = ActivationCoordinateActuator();
actu.setCoordinate(coord);
actu.setName('actuator_arm_flex_r');
actu.setOptimalForce(150);
actu.set_activation_time_constant(0.035);
actu.set_default_activation(0);
actu.set_min_control(-1);
actu.set_max_control(1);
model.addForce(actu);

coord = coords.get('arm_add_r');
actu = ActivationCoordinateActuator();
actu.setCoordinate(coord);
actu.setName('actuator_arm_add_r');
actu.setOptimalForce(150);
actu.set_activation_time_constant(0.035);
actu.set_default_activation(0);
actu.set_min_control(-1);
actu.set_max_control(1);
model.addForce(actu);

coord = coords.get('arm_rot_r');
actu = ActivationCoordinateActuator();
actu.setCoordinate(coord);
actu.setName('actuator_arm_rot_r');
actu.setOptimalForce(150);
actu.set_activation_time_constant(0.035);
actu.set_default_activation(0);
actu.set_min_control(-1);
actu.set_max_control(1);
model.addForce(actu);

% elbow left
coord = coords.get('elbow_flex_l');
actu = ActivationCoordinateActuator();
actu.setCoordinate(coord);
actu.setName('actuator_elbow_flex_l');
actu.setOptimalForce(150);
actu.set_activation_time_constant(0.035);
actu.set_default_activation(0);
actu.set_min_control(-1);
actu.set_max_control(1);
model.addForce(actu);

% elbow right
coord = coords.get('elbow_flex_r');
actu = ActivationCoordinateActuator();
actu.setCoordinate(coord);
actu.setName('actuator_elbow_flex_r');
actu.setOptimalForce(150);
actu.set_activation_time_constant(0.035);
actu.set_default_activation(0);
actu.set_min_control(-1);
actu.set_max_control(1);
model.addForce(actu);

model.finalizeConnections();
model.initSystem;
model.print([pathHere '\' osim_filename '_with_actuators.osim']);

%%

model = Model([pathHere '\' osim_filename '_with_actuators.osim']);
model.initSystem;

actuators = model.getActuators();
Nact = actuators.getSize();
actuator_info.coord_names = {};
actuator_info.coordi = [];
actuator_info.max_torque = [];
actuator_info.time_constant = [];

for i=1:Nact
    act_type = actuators.get(i-1).getConcreteClassName;
    if strcmp(act_type,'ActivationCoordinateActuator')
        act_i = actuators.get(i-1);
        acta_i = ActivationCoordinateActuator.safeDownCast(act_i);
        actuator_info.coord_names{end+1} = char(acta_i.getCoordinate().getName());
        actuator_info.max_torque(end+1) = double(acta_i.getOptimalForce());
        actuator_info.time_constant(end+1) = double(acta_i.get_activation_time_constant());
    end
end



















