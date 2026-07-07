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

if isempty(S.subject.IG_pelvis_y)
    import org.opensim.modeling.*;
    
    model = Model(osim_path);
    state = model.initSystem;
    
    calcn_or = model.getBodySet().get('calcn_r').findBaseFrame().getPositionInGround(state).getAsMat;
    pelvis_or = model.getBodySet().get('pelvis').findBaseFrame().getPositionInGround(state).getAsMat;
    
    model_info.IG_pelvis_y = pelvis_or(2)-calcn_or(2)+25e-3;

else
    model_info.IG_pelvis_y = S.subject.IG_pelvis_y;

end

end








