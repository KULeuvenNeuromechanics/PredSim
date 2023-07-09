function [contact_spheres] = get_contact_spheres(osim_path)
% --------------------------------------------------------------------------
% get_contact_spheres
%   This functions reads the parameter values describing the contact spheres 
%   of a given .osim model, and returns this information in a structured way.
% 
% INPUT:
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
% OUTPUT:
%   - contact_spheres -
%   * cell array of structs describing a contact sphere. 
%     Example cell:
%
%       % name of parent body
%       contact_spheres(1).body = 'calcn_r';
%       % name of contact sphere
%       contact_spheres(1).name = 's1_r';
%       % location in parent frame
%       contact_spheres(1).location = [0.0019 -0.01 -0.0038];
%       % radius of sphere
%       contact_spheres(1).radius = 0.032;
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
% Original date: 25/Nov/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


%% load model
import org.opensim.modeling.*;
model = Model(osim_path);
model.finalizeConnections();

forceSet = model.getForceSet();
geometrySet = model.getContactGeometrySet();


%%
csp = 1;

for i=0:forceSet.getSize()-1
    force_i = forceSet.get(i);
    if strcmp(force_i.getConcreteClassName(),"SmoothSphereHalfSpaceForce")

        force_i = SmoothSphereHalfSpaceForce.safeDownCast(force_i);

        socket1 = force_i.getSocket("sphere");
        socket1_objPath = char(socket1.getConnecteePath());
        [~,socket1_objName,~] = fileparts(socket1_objPath);
        geo = geometrySet.get(socket1_objName);
        geo1 = ContactSphere.safeDownCast(geo);
        geo1_body = geo1.getBody();

        % name of parent body
        contact_spheres(csp).body = char(geo1_body.getName());
        % name of contact sphere
        contact_spheres(csp).name = socket1_objName;
        % location in parent frame
        contact_spheres(csp).location = geo1.get_location().getAsMat();
        % radius of sphere
        contact_spheres(csp).radius = geo1.getRadius();
        % stiffness normal to plane
        contact_spheres(csp).stiffness = force_i.get_stiffness;
        % dissipation normal to plane
        contact_spheres(csp).dissipation = force_i.get_dissipation;
        % static friction coefficient in plane
        contact_spheres(csp).staticFriction = force_i.get_static_friction;
        % dynamic friction coefficient in plane
        contact_spheres(csp).dynamicFriction = force_i.get_dynamic_friction;
        % viscous friction coefficient in plane
        contact_spheres(csp).viscousFriction = force_i.get_viscous_friction;
        % transition velocity of static to dynamic friction
        contact_spheres(csp).transitionVelocity = force_i.get_transition_velocity;
        
        csp = csp+1;


    end




end













end