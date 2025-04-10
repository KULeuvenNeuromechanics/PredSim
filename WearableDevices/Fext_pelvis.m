function [ForceObj] = Fext_pelvis(init, settings_orthosis)

% create Orthosis object
ForceObj = Orthosis('ExtForce',init);

% get externalforce magnitude and point of applications
Fext = settings_orthosis.extForce;
r_origin = settings_orthosis.r_origin;

% try force proportional to pelvis velocity as an example
v_pevlvis_x = ForceObj.var_coord('pelvis_tx','vel');
%Fext_damper = -Fext * v_pevlvis_x; % example of a damper
Fext_damper = Fext * v_pevlvis_x; % example of something that injects energy

% add the body force
ForceObj.addBodyForce(Fext_damper, 'ext_force_pelvis',...
    'pelvis', r_origin, 'ground');

% include in outputs
ForceObj.addVarToPostProcessing(Fext_damper, 'Fpelvis')

end
