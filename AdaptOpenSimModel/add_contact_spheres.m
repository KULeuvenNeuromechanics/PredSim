function [] = add_contact_spheres(osim_path,contact_spheres)
% --------------------------------------------------------------------------
% add_contact_spheres
%   This functions adds contact spheres to the .osim model, and defines a
%   SmoothSphereHalfSpaceForce between each contact sphere and the ground.
% 
% INPUT:
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%   - contact_spheres -
%   * cell array of structs describing a contact sphere. 
%     Example cell:
%       % name of parent body
%       contact_spheres(1).body = 'calcn_r';
%       % name of contact sphere
%       contact_spheres(1).name = 's1_r';
%       % location in parent frame
%       contact_spheres(1).location = [0.0019 -0.01 -0.0038];
%       % radius of sphere
%       contact_spheres(1).radius = 0.032;
% 
%
% Original author: Lars D'Hondt
% Original date: 27/May/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


%% dynamic properties
stiffness           = 1000000;
dissipation         = 2.0;
staticFriction      = 0.8;
dynamicFriction     = 0.8;
viscousFriction     = 0.5;
transitionVelocity  = 0.2;

%% load model
import org.opensim.modeling.*;
model = Model(osim_path);


%% ground
ground_body = model.getGround();
groundContactSpace = ContactHalfSpace(Vec3(0,0,0),Vec3(0,0,-pi/2),ground_body);
groundContactSpace.setName('floor');
model.addContactGeometry(groundContactSpace);

%% loop over contact spheres
for i=1:length(contact_spheres)

% add contact geometry
    % get parent body
    body_i = model.getBodySet().get(contact_spheres(i).body);
    % construct sphere
    contact_sphere_1 = ContactSphere();
    % set radius
    contact_sphere_1.setRadius(contact_spheres(i).radius);
    % set location
    contactSphereLoc = Vec3.createFromMat(contact_spheres(i).location);
    contact_sphere_1.setLocation(contactSphereLoc);
    % assign to parent frame
    contact_sphere_1.setFrame(body_i);
    % set name
    contact_sphere_1.setName(contact_spheres(i).name);
    % add to model
    model.addContactGeometry(contact_sphere_1);   
    
% add contact force
    % construct force
    force_sphere_1 = SmoothSphereHalfSpaceForce();
    % set name
    force_sphere_1.setName(['SmoothSphereHalfSpaceForce_' contact_spheres(i).name]);
    % assign contact geometry and ground
    force_sphere_1.connectSocket_half_space(groundContactSpace)
    force_sphere_1.connectSocket_sphere(contact_sphere_1)
    % set dynamic properties
    force_sphere_1.set_stiffness(stiffness);
    force_sphere_1.set_dissipation(dissipation);
    force_sphere_1.set_static_friction(staticFriction);
    force_sphere_1.set_dynamic_friction(dynamicFriction);
    force_sphere_1.set_viscous_friction(viscousFriction);
    force_sphere_1.set_transition_velocity(transitionVelocity);
    % add to model
    model.addForce(force_sphere_1);
    
    % clear variable names
    clear 'contact_sphere_1' 'force_sphere_1'

end

%% save model
model.initSystem();
model.print(osim_path);






