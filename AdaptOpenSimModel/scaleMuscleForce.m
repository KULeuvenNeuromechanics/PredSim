function [] = scaleMuscleForce(osim_path_scaled, bodymass_generic, osim_path)
% --------------------------------------------------------------------------
% scaleMuscleForce
%   Scales the muscle properties of an OpenSim model that are not scaled
%   by the scale tool. This only includes properties that are used by
%   PredSim, i.e. 
%       max_isometric_force: ~mass^(2/3)
% 
%
% INPUT:
%   - osim_path_scaled -
%   * OpenSim model where muscles will be scaled. Muscles will be read
%   from this model, scaled, and written back to this model.
% 
%   - bodymass_generic -
%   * body mass of the generic (unscaled) model. Used to determine the mass
%   ratio for the scaling factor. Alternatively, pass the path to the
%   generic OpenSim model.
%
%   - osim_path - (optional) Default: osim_path_scaled
%   * OpenSim model before scaling muscles. Use this argument if you do not
%   want to overwrite the model file.
% 
%
% Original author: Lars D'Hondt
% Original date: 3 January 2025
% --------------------------------------------------------------------------

% set default for optional argument
if nargin <3
    osim_path = osim_path_scaled;
end

% get mass of generic model
if isnumeric(bodymass_generic) && numel(bodymass_generic)==1
    mass_generic = bodymass_generic;
else
    osim_path_generic = bodymass_generic;
    mass_generic = getModelMass(osim_path_generic);
end

% get mass of scaled model
mass_scaled = getModelMass(osim_path);

% scale factor
scale_factor = (mass_scaled/mass_generic)^(2/3);

import org.opensim.modeling.*;
model = Model(osim_path);

% apply scaling
for i=1:model.getMuscles().getSize()
    muscle_i = model.getMuscles().get(i-1);

    FMo_old = muscle_i.getMaxIsometricForce();
    FMo_new = FMo_old*scale_factor;
    muscle_i.setMaxIsometricForce(FMo_new);

end

% save model
model.initSystem();
model.print(osim_path_scaled);


end % end of function