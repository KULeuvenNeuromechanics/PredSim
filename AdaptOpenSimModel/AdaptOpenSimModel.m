% --------------------------------------------------------------------------
% AdaptOpenSimModel
%   This framework uses an OpenSim model file as input. Commonly used models 
%   do not contain all information needed to formulate the simulation. This
%   script adapts a given .osim file.
%
%   Settings:
%   - add_actuators_bool -
%   * Adds ideal torque actuators to the model.
%   See example below on how to specify actuator properties as a struct torq_act.
%
%   - add_contact_bool -
%   * Adds contact spheres that interact with the ground.
%   See example below on how to specify contact spheres properties as a
%   struct contact_spheres.
%
%   - use_reference_contacts_bool -
%   * Creates a struct contact_spheres based on the contact properties of a
%   reference model.
%
%   - scale_contact_location_bool -
%   * OpenSim scale tool does not adapt contact sphere locations. Use this
%   setting to scale contact sphere locations with respect to a reference
%   model.
% 
% Original author: Lars D'Hondt
% Original date: 27/May/2022
%
% Last edit by: Bram Van Den Bosch  
% Last edit date: 25/January/2024
% --------------------------------------------------------------------------


clear
clc
[pathHere,~,~] = fileparts(mfilename('fullpath'));

% .osim file to adapt
path_osim_in = fullfile(pathHere,'Falisse_et_al_2022.osim');

% adapted .osim file
path_osim_out = fullfile(pathHere,'Falisse_et_al_2022.osim');

% reference model for contact
path_reference_model = fullfile(pathHere,'Falisse_et_al_2022.osim');

% select what to do
add_actuators_bool = 1;
add_contact_bool = 0;
use_reference_contacts_bool = 1;
scale_contact_spheres = 1;
scale_contact_location_bool = 1;

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

%% Define scaling factors for contact spheres

scale.stiffness = 1;
scale.dissipation = 1;

%%

import org.opensim.modeling.*;
model = Model(path_osim_in);
model_name = char(model.getName);
model_name = [model_name '_AdaptedForPredSim'];
model.setName(model_name);
model.print(path_osim_out);


%%
if use_reference_contacts_bool
    contact_spheres = get_contact_spheres(path_reference_model);
end
if add_actuators_bool
    add_actuators(path_osim_out,torq_act);
end
if add_contact_bool
    add_contact_spheres(path_osim_out,contact_spheres);
end
if scale_contact_spheres
    scaleContactSpheres(path_reference_model,path_osim_out,path_osim_out,scale)
end
if scale_contact_location_bool
    fixContactSpherePositionAfterScaling(path_reference_model,path_osim_out);
end

