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
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


coordinates = fields(model_info.ExtFunIO.coordi);
for i = 1:length(model_info.ExtFunIO.jointi.(group))
    for c = 1:length(coordinates)
        if model_info.ExtFunIO.jointi.(group)(i) == model_info.ExtFunIO.coordi.(coordinates{c})
            model_info.ExtFunIO.coord_names.(group){i} = coordinates{c};
        end
    end
end
