function [] = add_contact_spheres(osim_path,contact_spheres,varargin)
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
%
%       Minimal input fields - contact geometry parameters
%       % name of parent body
%       contact_spheres(1).body = 'calcn_r';
%       % name of contact sphere
%       contact_spheres(1).name = 's1_r';
%       % location in parent frame
%       contact_spheres(1).location = [0.0019 -0.01 -0.0038];
%       % radius of sphere
%       contact_spheres(1).radius = 0.032;
%
%       Optional input fields - mechanical properties
%       % stiffness normal to plane
%       contact_spheres(1).stiffness = 1000000;
%       % dissipation normal to plane
%       contact_spheres(1).dissipation = 2.0;
%       % static friction coefficient in plane
%       contact_spheres(1).staticFriction = 0.8;
%       % dynamic friction coefficient in plane
%       contact_spheres(1).dynamicFriction = 0.8;
%       % viscous friction coefficient in plane
%       contact_spheres(1).viscousFriction = 0.5;
%       % transition velocity of static to dynamic friction
%       contact_spheres(1).transitionVelocity = 0.2;
%       
% 
%
% Original author: Lars D'Hondt
% Original date: 27/May/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


%% dynamic properties
default_stiffness           = 1000000;
default_dissipation         = 2.0;
default_staticFriction      = 0.8;
default_dynamicFriction     = 0.8;
default_viscousFriction     = 0.5;
default_transitionVelocity  = 0.2;

%% extract optional inputs
for i=1:2:length(varargin)
    if strcmp(varargin{i},stiffness)
        stiffness = varargin{i+1};
    end
end


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
    % Look for optional inputs of dynamic properties
    if isfield(contact_spheres(i),'stiffness')
        stiffness = contact_spheres(i).stiffness;
    else
        stiffness = default_stiffness;
    end
    if isfield(contact_spheres(i),'dissipation')
        dissipation = contact_spheres(i).dissipation;
    else
        dissipation = default_dissipation;
    end
    if isfield(contact_spheres(i),'staticFriction')
        staticFriction = contact_spheres(i).staticFriction;
    else
        staticFriction = default_staticFriction;
    end
    if isfield(contact_spheres(i),'dynamicFriction')
        dynamicFriction = contact_spheres(i).dynamicFriction;
    else
        dynamicFriction = default_dynamicFriction;
    end
    if isfield(contact_spheres(i),'viscousFriction')
        viscousFriction = contact_spheres(i).viscousFriction;
    else
        viscousFriction = default_viscousFriction;
    end
    if isfield(contact_spheres(i),'transitionVelocity')
        transitionVelocity = contact_spheres(i).transitionVelocity;
    else
        transitionVelocity = default_transitionVelocity;
    end

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
model.finalizeConnections();
model.initSystem();
model.print(osim_path);






