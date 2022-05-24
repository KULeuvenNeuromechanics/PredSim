

% Lars D'Hondt
% 24/May/2022


clear
clc


% contact spheres right side
csp = 1;

contact_spheres(csp).body = 'calcn_r';
contact_spheres(csp).name = 's1_r';
contact_spheres(csp).location = [0.0019 -0.01 -0.0038];
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




%%


import org.opensim.modeling.*;
[pathHere,~,~] = fileparts(mfilename('fullpath'));

osim_filename = 'Rajagopal2015_opensense';

model = Model([pathHere '\' osim_filename '.osim']);


% ground
ground_body = model.getGround();
groundContactSpace = ContactHalfSpace(Vec3(0,0,0),Vec3(0,0,-1.57),ground_body);
groundContactSpace.setName('floor');
model.addContactGeometry(groundContactSpace);

for i=1:length(contact_spheres)

    % spheres
    body_i = model.getBodySet().get(contact_spheres(i).body);
    contactSphereRadius = contact_spheres(i).radius;
    contactSphereLoc = Vec3.createFromMat(contact_spheres(i).location);
    
    contact_sphere_1 = ContactSphere();
    contact_sphere_1.setRadius(contactSphereRadius);
    contact_sphere_1.setLocation(contactSphereLoc);
    contact_sphere_1.setFrame(body_i);
    contact_sphere_1.setName(contact_spheres(i).name);
    model.addContactGeometry(contact_sphere_1);
    
    
    %%
    stiffness           = 1000000;
    dissipation         = 2.0;
    staticFriction      = 0.8;
    dynamicFriction     = 0.8;
    viscousFriction     = 0.5;
    transitionVelocity  = 0.2;
    
    force_sphere_1 = SmoothSphereHalfSpaceForce();
    force_sphere_1.setName(['SmoothSphereHalfSpaceForce_' contact_spheres(i).name]);
    force_sphere_1.connectSocket_half_space(groundContactSpace)
    force_sphere_1.connectSocket_sphere(contact_sphere_1)
    force_sphere_1.set_stiffness(stiffness);
    force_sphere_1.set_dissipation(dissipation);
    force_sphere_1.set_static_friction(staticFriction);
    force_sphere_1.set_dynamic_friction(dynamicFriction);
    force_sphere_1.set_viscous_friction(viscousFriction);
    force_sphere_1.set_transition_velocity(transitionVelocity);
    model.addForce(force_sphere_1);
    
    clear 'contact_sphere_1' 'force_sphere_1'


end
%%
model.initSystem();

model.print([pathHere '\' osim_filename '_test.osim']);






