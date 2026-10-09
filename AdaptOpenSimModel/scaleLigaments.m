function [] = scaleLigaments(osim_path_scaled, bodymass_generic)
% --------------------------------------------------------------------------
% scaleLigaments
%   Scales the ligament properties of an OpenSim model that are not scaled
%   by the scale tool. This only includes properties that are used by
%   PredSim, i.e. 
%       pcsa_force: ~mass^(2/3)
% 
%
% INPUT:
%   - osim_path_scaled -
%   * OpenSim model where ligaments will be scaled. Ligaments will be read
%   from this model, scaled, and written back to this model.
% 
%   - bodymass_generic -
%   * body mass of the generic (unscaled) model. Used to determine the mass
%   ratio for the scaling factor. Alternatively, pass the path to the
%   generic OpenSim model.
%
% 
% Original author: Lars D'Hondt
% Original date: 9 February 2024
%
% --------------------------------------------------------------------------
% This file is part of PredSim.
% 
% PredSim: A Framework for Rapid Predictive Simulations of Locomotion
% Copyright (c) 2026 KU Leuven
% 
% PredSim is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Affero General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% 
% PredSim is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public 
% License for more details.
% 
% You should have received a copy of the GNU Affero General Public License 
% along with PredSim. If not, see <https://www.gnu.org/licenses/>.
% --------------------------------------------------------------------------

% get mass of generic model
if isnumeric(bodymass_generic) && numel(bodymass_generic)==1
    mass_generic = bodymass_generic;
else
    osim_path_generic = bodymass_generic;
    mass_generic = getModelMass(osim_path_generic);
end

% get mass of scaled model
mass_scaled = getModelMass(osim_path_scaled);

% scale factor
scale_factor = (mass_scaled/mass_generic)^(2/3);

import org.opensim.modeling.*;
model = Model(osim_path_scaled);

% apply scaling
for i=1:model.getForceSet().getSize()
    force_i = model.getForceSet().get(i-1);
    if strcmp(force_i.getConcreteClassName(),'Ligament')
        lig = Ligament.safeDownCast(force_i);
        PCSA_force_old = lig.get_pcsa_force();
        PCSA_force_new = PCSA_force_old *scale_factor;
        lig.set_pcsa_force(PCSA_force_new);

    end
end

% save model
model.initSystem();
model.print(osim_path_scaled);


end % end of function