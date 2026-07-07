function model_info = addCoordNames(model_info,group)
% --------------------------------------------------------------------------
% addCoordNames
%   This function creates a field with the coordinate names corresponding
%   to a subset of indices.
%   
% INPUT:
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
%   - group -
%   * subset of indices
%
% OUTPUT:
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% Original author: Dhruv Gupta
% Original date: 
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


coordinates = fields(model_info.ExtFunIO.coordi);
for i = 1:length(model_info.ExtFunIO.jointi.(group))
    for c = 1:length(coordinates)
        if model_info.ExtFunIO.jointi.(group)(i) == model_info.ExtFunIO.coordi.(coordinates{c})
            model_info.ExtFunIO.coord_names.(group){i} = coordinates{c};
        end
    end
end
