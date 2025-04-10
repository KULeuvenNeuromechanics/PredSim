function [ForceObj] = Fext_pelvis_emg(init, settings_orthosis)

% create Orthosis object
ForceObj = Orthosis('ExtForce',init);

% get externalforce magnitude and point of applications
Fext = settings_orthosis.extForce;
r_origin = settings_orthosis.r_origin;

% force as a function of the gait cycle
emg_sum = 0;
for i = 1:length(settings_orthosis.muscle_names)
    emg_sum = emg_sum + ForceObj.var_muscle(settings_orthosis.muscle_names{i});
end
emg_mean = emg_sum./length(settings_orthosis.muscle_names);

F_pelvis = emg_mean*Fext;
% add the body force
ForceObj.addBodyForce(F_pelvis, 'ext_force_pelvis',...
    'pelvis', r_origin, 'ground');

% include in outputs
ForceObj.addVarToPostProcessing(F_pelvis, 'Fpelvis')

end
