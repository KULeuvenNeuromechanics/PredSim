function [scaling] = getScaleFactor(S,model_info,bounds_nsc)
% --------------------------------------------------------------------------
% getScaleFactor
%   This script provides scaling factors for the optimisation variables.
%   Scale factors are based on the bound with highest absolute value, so
%   optimisation variables remain within the interval [-1,1].
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
%   - bounds_nsc -
%   * boundaries for all optimisation variables (not scaled)
%
% OUTPUT:
%   - scaling -
%   * scale factors for all optimisation variables
% 
% Original author: Lars D'Hondt
% Original date: 6/April/2023
%
% Last edit by:
% Last edit date:
% --------------------------------------------------------------------------


coordinate_names = model_info.ExtFunIO.coord_names.all;
NCoord = model_info.ExtFunIO.jointi.nq.all;

%% Default: based on bounds
% Qs
scaling.Qs = max(abs(bounds_nsc.Qs.lower),abs(bounds_nsc.Qs.upper));
% Qdots
scaling.Qdots = max(abs(bounds_nsc.Qdots.lower),abs(bounds_nsc.Qdots.upper));
% Qdotdots
scaling.Qdotdots = max(abs(bounds_nsc.Qdotdots.lower),abs(bounds_nsc.Qdotdots.upper));
% Torque actuator 
scaling.ActuatorTorque = struct_array_to_double_array(model_info.actuator_info.parameters,'max_torque');
% Time derivative of muscle activations
scaling.vA = 100;
% Muscle activations
scaling.a = 1;
% Torque actuator activations
scaling.a_a = 1;
% Torque actuator excitations
scaling.e_a = 1;
% Time derivative of muscle-tendon forces
scaling.dFTtilde = 100;
% Muscle-tendon forces
scaling.FTtilde = max(abs(bounds_nsc.FTtilde.lower),abs(bounds_nsc.FTtilde.upper)); 
% Joint moment equality constraints
scaling.Moments = ones(1,NCoord);
% Kinematic constraints contact forces
scaling.F_constr = 100 * model_info.mass*9.81;

%% Based on provided table
fields = ["Qs","Qdots","Qdotdots","Moments"];
if exist(S.misc.default_scaling_NLP,'file')
    default_scaling = readtable(S.misc.default_scaling_NLP);
    for j=1:NCoord
        default_scaling_j = default_scaling(strcmp(default_scaling.name,coordinate_names{j}),:);
        if ~isempty(default_scaling_j)
            for i=fields
                scaling.(i)(j) = default_scaling_j.(char(i));    
            end
        end
    end
end

%% Adjust scaling based on settings
for i=fields
    if ~isempty(S.misc.(['scaling_' char(i)]))
        [new_sf] = unpack_name_value_combinations(S.misc.(['scaling_' char(i)]),coordinate_names,1);

        for j=1:NCoord
            coordinate = coordinate_names{j};
            coord_idx = model_info.ExtFunIO.coordi.(coordinate);
            if ~isnan(new_sf(j))
                scaling.(i)(coord_idx) = new_sf(j);
            end
        end
    end
end


end % end of function