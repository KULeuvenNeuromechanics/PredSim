function [ForceObj] = Fext_torso_bumpm(init, settings_orthosis)

% create Orthosis object
ForceObj = Orthosis('ExtForce',init);

% get externalforce magnitude and point of applications
Fext = settings_orthosis.extForce;
r_origin = settings_orthosis.r_origin;

% add the body force
ForceObj.addBodyForce(Fext, 'ext_force_pelvis',...
    'torso', r_origin, 'ground');

% include in outputs
ForceObj.addVarToPostProcessing(Fext, 'Ftorso')

end
