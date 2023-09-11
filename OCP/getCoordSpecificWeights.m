function [S] = getCoordSpecificWeights(S,model_info)
% --------------------------------------------------------------------------
% getCoordSpecificWeights
%   This script provides coordinate specific weights, such that the
%   contribution of each coordinate acceleration and passive torque to the
%   relevant cost function terms can be selected.
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - S -
%   * setting structure S
% 
% Original author: Lars D'Hondt
% Original date: 11/September/2023
%
% Last edit by:
% Last edit date:
% --------------------------------------------------------------------------

coordinate_names = model_info.ExtFunIO.coord_names.all;
NCoord = model_info.ExtFunIO.jointi.nq.all;

%% Initialise with 1
S.weights.coordspecific_Qdotdots = ones(1,NCoord);
S.weights.coordspecific_PassTorq = ones(1,NCoord);

%% Based on provided table
fields = ["Qdotdots","PassTorq"];
if exist(S.weights.default_coordspecific_weights,'file')
    default_weight = readtable(S.weights.default_coordspecific_weights);
    for j=1:NCoord
        default_weight_j = default_weight(strcmp(default_weight.name,coordinate_names{j}),:);
        if ~isempty(default_weight_j)
            for i=fields
                S.weights.(['coordspecific_' char(i)])(j) = default_weight_j.(char(i));    
            end
        end
    end
end

%% Adjust scaling based on settings
for i=fields
    if ~isempty(S.weights.(['set_coordspecific_' char(i)]))
        [new_w] = unpack_name_value_combinations(S.weights.(['set_coordspecific_' char(i)]),coordinate_names,1);

        for j=1:NCoord
            coordinate = coordinate_names{j};
            coord_idx = model_info.ExtFunIO.coordi.(coordinate);
            if ~isnan(new_w(j))
                S.weights.(['coordspecific_' char(i)])(coord_idx) = new_w(j);
            end
        end
    end
end

% exclude PassTorq weights for coordinates that have no limit torque
S.weights.coordspecific_PassTorq = S.weights.coordspecific_PassTorq(model_info.ExtFunIO.jointi.limitTorque);

end % end of function