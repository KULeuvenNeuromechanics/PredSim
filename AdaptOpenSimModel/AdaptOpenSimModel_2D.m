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
%   - scale_contact_spheres_bool -
%   * Scales the contact spheres' radii based on foot length
%   * Scales the contact spheres' stiffness and dissipation based on provided struct (see example below) 
%
%   - scale_contact_location_bool -
%   * OpenSim scale tool does not adapt contact sphere locations. Use this
%   setting to scale contact sphere locations with respect to a reference
%   model.
% 
%   - scale_ligaments_bool -
%   * OpenSim scale tool does not adapt ligament pcsa. Use this setting to 
%   scale ligament pcsa_force with respect to a reference model.
%
%
% Original author: Lars D'Hondt
% Original date: 27/May/2022
%
% Last edit by: Ellis van Can  
% Last edit date: 12/November/2025
% --------------------------------------------------------------------------

clear
close all;
clc
% [pathHere,~,~] = fileparts(mfilename('fullpath'));

% .osim file to adapt
SubjectName = 'test_bug_2D';
path_osim = "C:\GBW_MyPrograms\PredSimSHARED\Tests\BugContactSpheres";
path_osim_in = fullfile(path_osim,[SubjectName,'.osim']);

% adapted .osim file
path_osim_out = fullfile(path_osim,[SubjectName,'_Adapted.osim']);

% reference model for contact
path_reference_model = 'C:\GBW_MyPrograms\PredSim\Subjects\gait1018\gait1018.osim';

% select what to do
scaleFMO = 0;
add_actuators_bool = 0;
add_contact_bool = 0;
use_reference_contacts_bool = 0;
scale_contact_spheres_bool = 1;
scale_contact_location_bool = 1;
scale_ligaments_bool = 1;

scale.stiffness = 0.8;
scale.dissipation = 0.98;

factorR = 0.71;
factorL = 0.69;

mass_subject = 42.3;
mass_genmodel = 62;
%% Define contact spheres

% contact spheres right side
csp = 1;

% name of parent body
contact_spheres(csp).body = 'calcn_r';
% name of contact sphere
contact_spheres(csp).name = 'heel_r';
% location in parent frame
contact_spheres(csp).location = [0.031307527581931796 0.010435842527310599 0];
% radius of sphere
contact_spheres(csp).radius = 0.035;
csp = csp+1;

contact_spheres(csp).body = 'calcn_r';
contact_spheres(csp).name = 'front_r';
contact_spheres(csp).location = [0.17740932296428019 -0.015653763790965898 0.0052179212636552993];
contact_spheres(csp).radius = 0.015;
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


% %% Define ideal torque actuators
% 
% ita = 1;
% 
% % coordinate to which te actuator applies torque
% torq_act(ita).coord = ('arm_flex_r');
% % maximum torque
% torq_act(ita).max_torque = 150;
% % time constant of the activation dynamics
% torq_act(ita).time_constant = 0.035;
% ita = ita+1;
% 
% % coordinate to which te actuator applies torque
% torq_act(ita).coord = ('arm_add_r');
% % maximum torque
% torq_act(ita).max_torque = 150;
% % time constant of the activation dynamics
% torq_act(ita).time_constant = 0.035;
% ita = ita+1;
% 
% % coordinate to which te actuator applies torque
% torq_act(ita).coord = ('arm_rot_r');
% % maximum torque
% torq_act(ita).max_torque = 150;
% % time constant of the activation dynamics
% torq_act(ita).time_constant = 0.035;
% ita = ita+1;
% 
% % coordinate to which te actuator applies torque
% torq_act(ita).coord = ('elbow_flex_r');
% % maximum torque
% torq_act(ita).max_torque = 150;
% % time constant of the activation dynamics
% torq_act(ita).time_constant = 0.035;
% ita = ita+1;
% 
% % mirror to get left side
% for i=1:length(torq_act)
%     torq_act(ita).coord = [torq_act(i).coord(1:end-1) 'l'];
%     torq_act(ita).max_torque = torq_act(i).max_torque;
%     torq_act(ita).time_constant = torq_act(i).time_constant;
%     ita = ita+1;
% end


%%

import org.opensim.modeling.*;
model = Model(path_osim_in);
model_name = char(model.getName);
model_name = [model_name '_AdaptedForPredSim'];
model.setName(model_name);
model.print(path_osim_out);


%%
if scaleFMO
    if ~contains(model_name,'_sf')
    mass_genmodel = 62;
    adjust_FMo_and_print(path_osim_in,mass_genmodel,mass_subject,path_osim_out);
    end
end

SimmSpline_joint_to_polynomial(path_osim_out);

if use_reference_contacts_bool
    contact_spheres = get_contact_spheres(path_reference_model);
end
if add_actuators_bool
    add_actuators(path_osim_out,torq_act);
end
if add_contact_bool
    add_contact_spheres(path_osim_out,contact_spheres);
end
if scale_contact_spheres_bool
    scaleContactSpheres(path_reference_model,path_osim_in,path_osim_out,scale,factorR,factorL)
end
if scale_contact_location_bool
    fixContactSpherePositionAfterScaling(path_reference_model,path_osim_out,factorR,factorL);
end
if scale_ligaments_bool
    scaleLigaments(path_osim_out, path_reference_model);
end

