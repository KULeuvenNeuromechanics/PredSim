% --------------------------------------------------------------------------
% AdaptOpenSimModel
%   This framework uses an OpenSim model file as input. Commonly used models 
%   do not contain all information needed to formulate the simulation. This
%   script adapts a given .osim file to include:
%   1) Ideal torque actuators
%   2) Contact elements
% 
% 
% Original author: Lars D'Hondt
% Original date: 27/May/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


clear
clc
[pathHere,~,~] = fileparts(mfilename('fullpath'));

% .osim file to adapt
path_osim_in = fullfile(pathHere,'Falisse_et_al_2022.osim');

% adapted .osim file
path_osim_out = fullfile(pathHere,'Falisse_et_al_2022.osim');

add_actuators_bool = 1;
add_contact_bool = 0;

%% Define contact spheres

% contact spheres right side
csp = 1;

% name of parent body
contact_spheres(csp).body = 'calcn_r';
% name of contact sphere
contact_spheres(csp).name = 's1_r';
% location in parent frame
contact_spheres(csp).location = [0.0019 -0.01 -0.0038];
% radius of sphere
contact_spheres(csp).radius = 0.032;
csp = csp+1;

contact_spheres(csp).body = 'calcn_r';
contact_spheres(csp).name = 's2_r';
contact_spheres(csp).location = [0.1483 -0.01 -0.0287];
contact_spheres(csp).radius = 0.032;
csp = csp+1;

contact_spheres(csp).body = 'calcn_r';
contact_spheres(csp).name = 's3_r';
contact_spheres(csp).location = [0.1330 -0.01 0.0516];
contact_spheres(csp).radius = 0.032;
csp = csp+1;

contact_spheres(csp).body = 'calcn_r';
contact_spheres(csp).name = 's4_r';
contact_spheres(csp).location = [0.06623 -0.01 0.02636];
contact_spheres(csp).radius = 0.032;
csp = csp+1;

contact_spheres(csp).body = 'toes_r';
contact_spheres(csp).name = 's5_r';
contact_spheres(csp).location = [0.06 -0.01 -0.01876];
contact_spheres(csp).radius = 0.032;
csp = csp+1;

contact_spheres(csp).body = 'toes_r';
contact_spheres(csp).name = 's6_r';
contact_spheres(csp).location = [0.045 -0.01 0.06186];
contact_spheres(csp).radius = 0.032;
csp = csp+1;

% mirror to get left side
for i=1:length(contact_spheres)
    contact_spheres(csp).body = [contact_spheres(i).body(1:end-1) 'l'];
    contact_spheres(csp).name = [contact_spheres(i).name(1:end-1) 'l'];
    contact_spheres(csp).location = contact_spheres(i).location;
    contact_spheres(csp).location(3) = -contact_spheres(csp).location(3);
    contact_spheres(csp).radius = contact_spheres(i).radius;
    csp = csp+1;
end


%% Define ideal torque actuators

ita = 1;

% coordinate to which te actuator applies torque
torq_act(ita).coord = ('arm_flex_r');
% maximum torque
torq_act(ita).max_torque = 150;
% time constant of the activation dynamics
torq_act(ita).time_constant = 0.035;
ita = ita+1;

% coordinate to which te actuator applies torque
torq_act(ita).coord = ('arm_add_r');
% maximum torque
torq_act(ita).max_torque = 150;
% time constant of the activation dynamics
torq_act(ita).time_constant = 0.035;
ita = ita+1;

% coordinate to which te actuator applies torque
torq_act(ita).coord = ('arm_rot_r');
% maximum torque
torq_act(ita).max_torque = 150;
% time constant of the activation dynamics
torq_act(ita).time_constant = 0.035;
ita = ita+1;

% coordinate to which te actuator applies torque
torq_act(ita).coord = ('elbow_flex_r');
% maximum torque
torq_act(ita).max_torque = 150;
% time constant of the activation dynamics
torq_act(ita).time_constant = 0.035;
ita = ita+1;

% mirror to get left side
for i=1:length(torq_act)
    torq_act(ita).coord = [torq_act(i).coord(1:end-1) 'l'];
    torq_act(ita).max_torque = torq_act(i).max_torque;
    torq_act(ita).time_constant = torq_act(i).time_constant;
    ita = ita+1;
end


%%

import org.opensim.modeling.*;
model = Model(path_osim_in);
model_name = char(model.getName);
model_name = [model_name '_AdaptedForPredSim'];
model.setName(model_name);
model.print(path_osim_out);


%%
if add_actuators_bool
    add_actuators(path_osim_out,torq_act);
end
if add_contact_bool
    add_contact_spheres(path_osim_out,contact_spheres);
    
%     Fix the location of contact spheres based on the size of model
    fixContactSpherePositionAfterScaling('../Subjects/Falisse_et_al_2022/Falisse_et_al_2022.osim',path_osim_out);
end

