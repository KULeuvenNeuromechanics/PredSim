function [model_info] = get_IG_pelvis_y(S,osim_path,model_info)
% --------------------------------------------------------------------------
% get_IG_pelvis_y
%   This function returns the pelvis height for quasi-random initial guess
%
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%       
% 
% Original author: Lars D'Hondt
% Original date: 12/April/2022
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

import org.opensim.modeling.*;

name_floating_base_vertical = model_info.ExtFunIO.coord_names.all{...
        model_info.ExtFunIO.jointi.base_vertical};

if ~isempty(S.subject.IG_pelvis_y) && isnumeric(S.subject.IG_pelvis_y)
    
    model_info.IG_pelvis_y = S.subject.IG_pelvis_y;

    value_source = 'S.subject.IG_pelvis_y';

elseif ~isempty(S.subject.IG_pelvis_y) ...
        && strcmpi(S.subject.IG_pelvis_y,'default_pose')

    model = Model(osim_path);
    
    coord_floating_base_vertical = model.getCoordinateSet().get(...
        name_floating_base_vertical);
    model_info.IG_pelvis_y = coord_floating_base_vertical.getDefaultValue;

    value_source = 'OpenSim model: default pose';

else
    model = Model(osim_path);
    state = model.initSystem;

    try
        calcn_or = model.getBodySet().get('calcn_r').findBaseFrame().getPositionInGround(state).getAsMat;
        pelvis_or = model.getBodySet().get('pelvis').findBaseFrame().getPositionInGround(state).getAsMat;
    catch
        error("Unable to automatically determine initial guess for %s.\n" + ...
            "Please use S.subject.IG_pelvis_y to set the value, or " + ...
            "use S.subject.IG_pelvis_y = 'default_pose' to use the default " + ...
            "pose value from the OpenSim model.")
    end

    model_info.IG_pelvis_y = pelvis_or(2)-calcn_or(2)+25e-3;

    value_source = 'OpenSim model: calcn_r-pelvis y-distance';
  
end

% Inform user which value is used and where it comes from
if strcmp(S.solver.IG_selection, 'quasi-random')
    fprintf("\tUsing %5.3fm as initial guess for %s. (%s)\n",...
        model_info.IG_pelvis_y, name_floating_base_vertical, value_source);

elseif S.subject.adapt_IG_pelvis_y
    fprintf("\tAdapting initial guess to use %5.3fm as mean value for %s. (%s)\n",...
        model_info.IG_pelvis_y, name_floating_base_vertical, value_source);
end

end








