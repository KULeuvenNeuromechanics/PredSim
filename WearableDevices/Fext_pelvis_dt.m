function [ForceObj] = Fext_pelvis_dt(init, settings_orthosis)

% create Orthosis object
ForceObj = Orthosis('ExtForce',init, true);

% get externalforce magnitude and point of applications
Fext = settings_orthosis.extForce;
r_origin = settings_orthosis.r_origin;

% force as a function of the gait cycle
N_control = ForceObj.getNmesh();

dt_vect = linspace(0,1,N_control);
F_sin = Fext * sin(dt_vect*2*pi);

figure();
plot(F_sin');

% add the body force
ForceObj.addBodyForce(F_sin, 'ext_force_pelvis',...
    'pelvis', r_origin, 'ground');

% include in outputs
ForceObj.addVarToPostProcessing(F_sin, 'Fpelvis')

end
