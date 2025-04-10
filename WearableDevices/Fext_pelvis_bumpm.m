function [ForceObj] = Fext_pelvis_bumpm(init, settings_orthosis)

% create Orthosis object
ForceObj = Orthosis('ExtForce',init);

% get externalforce magnitude and point of applications
Fext = settings_orthosis.extForce;
r_origin = settings_orthosis.r_origin;

% add the body force
ForceObj.addBodyForce(Fext, 'ext_force_pelvis',...
    'pelvis', r_origin, 'ground');

% include in outputs
ForceObj.addVarToPostProcessing(Fext, 'Fpelvis')

end
