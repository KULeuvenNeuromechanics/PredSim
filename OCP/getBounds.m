function [bounds_nsc] = getBounds(S,model_info)
% --------------------------------------------------------------------------
% getBounds
%   This script provides bounds for the (not scaled) optimisation variables.
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - bounds_nsc -
%   * boundaries for all optimisation variables (not scaled)
%
% 
% Original author: Lars D'Hondt
% Original date: 5/April/2023
%
% Last edit by:
% Last edit date:
% --------------------------------------------------------------------------


% Get the names of the coordinates
coordinate_names = model_info.ExtFunIO.coord_names.all;
NCoord = model_info.ExtFunIO.jointi.nq.all;
NMuscle = model_info.muscle_info.NMuscle;

%% Initialise
coords = ["Qs","Qdots","Qdotdots"];
for i=coords
    bounds_nsc.(i).lower = nan(1,NCoord);
    bounds_nsc.(i).upper = nan(1,NCoord);
end

%% Get default bounds
if exist(S.bounds.default_coordinate_bounds,'file')
    default_bounds = readtable(S.bounds.default_coordinate_bounds);
    for j=1:NCoord
        default_bounds_j = default_bounds(strcmp(default_bounds.name,coordinate_names{j}),:);
        if ~isempty(default_bounds_j)
            for i=coords
                bounds_nsc.(i).lower(j) = default_bounds_j.([char(i) '_lower']);
                bounds_nsc.(i).upper(j) = default_bounds_j.([char(i) '_upper']);    
            end
        end
    end
end

%% Read bounds from osim file
import org.opensim.modeling.*;
model = Model(model_info.osim_path);

for j=1:NCoord
    coord_j = model.getCoordinateSet().get(coordinate_names{j});
    lb_j = coord_j.getRangeMin();
    ub_j = coord_j.getRangeMax();
    
    % Qs
    if isnan(bounds_nsc.Qs.lower(j))
        bounds_nsc.Qs.lower(j) = lb_j;
    end
    if isnan(bounds_nsc.Qs.upper(j))
        bounds_nsc.Qs.upper(j) = ub_j;
    end
    % Range of motion
    range_j = bounds_nsc.Qs.upper(j) - bounds_nsc.Qs.lower(j);
    % Qdots
    if isnan(bounds_nsc.Qdots.lower(j))
        bounds_nsc.Qdots.lower(j) = -range_j*S.bounds.Qdots_factor_RoM;
    end
    if isnan(bounds_nsc.Qdots.upper(j))
        bounds_nsc.Qdots.upper(j) = range_j*S.bounds.Qdots_factor_RoM;
    end
    % Qdotdots
    if isnan(bounds_nsc.Qdotdots.lower(j))
        bounds_nsc.Qdotdots.lower(j) = -range_j*S.bounds.Qdotdots_factor_RoM;
    end
    if isnan(bounds_nsc.Qdotdots.upper(j))
        bounds_nsc.Qdotdots.upper(j) = range_j*S.bounds.Qdotdots_factor_RoM;
    end


end


%% Adjust bounds based on settings
for i=coords
    if ~isempty(S.bounds.(i))
        [new_lb,new_ub] = unpack_name_value_combinations(S.bounds.(i),coordinate_names,[1,1]);
        
        new_lb(model_info.ExtFunIO.jointi.rotations) = new_lb(model_info.ExtFunIO.jointi.rotations)*pi/180;
        new_ub(model_info.ExtFunIO.jointi.rotations) = new_ub(model_info.ExtFunIO.jointi.rotations)*pi/180;
        
        for j=1:NCoord
            coordinate = coordinate_names{j};
            coord_idx = model_info.ExtFunIO.coordi.(coordinate);
        
            if ~isnan(new_lb(j))
                bounds_nsc.(i).lower(coord_idx) = new_lb(j);
            end
        
            if ~isnan(new_ub(j))
                bounds_nsc.(i).upper(coord_idx) = new_ub(j);
            end
        end
    
    end
end


%% Vertical position of floating base
if ~isempty(S.bounds.factor_IG_pelvis_ty.lower)
    bounds_nsc.Qs.lower(model_info.ExtFunIO.jointi.base_vertical) = model_info.IG_pelvis_y...
        *S.bounds.factor_IG_pelvis_ty.lower;
end
if ~isempty(S.bounds.factor_IG_pelvis_ty.upper)
    bounds_nsc.Qs.upper(model_info.ExtFunIO.jointi.base_vertical) = model_info.IG_pelvis_y...
        *S.bounds.factor_IG_pelvis_ty.upper;
end


%% Adjust for half gait cycle
if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    bounds_nsc.Qs.upper(model_info.ExtFunIO.jointi.base_forward) = ...
        bounds_nsc.Qs.upper(model_info.ExtFunIO.jointi.base_forward)/2;
end

%% Hard bounds
% We impose the initial position of pelvis_tx to be 0
bounds_nsc.Qs_0.lower = bounds_nsc.Qs.lower;
bounds_nsc.Qs_0.upper = bounds_nsc.Qs.upper;
bounds_nsc.Qs_0.lower(model_info.ExtFunIO.jointi.base_forward) = 0;
bounds_nsc.Qs_0.upper(model_info.ExtFunIO.jointi.base_forward) = 0;


%% Muscle activations
bounds_nsc.a.lower = S.bounds.activation_all_muscles.lower*ones(1,NMuscle);
bounds_nsc.a.upper = ones(1,NMuscle);

if ~isempty(S.bounds.activation_selected_muscles)
    [new_lb,new_ub] = unpack_name_value_combinations(S.bounds.activation_selected_muscles,...
        muscle_info.muscle_names,[1,1]);
    for j=1:NMuscle
        if ~isnan(new_lb(j))
            bounds_nsc.a.lower(j) = new_lb(j);
        end
        if ~isnan(new_ub(j))
            bounds_nsc.a.upper(j) = new_ub(j);
        end
    end 
end

%% Final time
bounds_nsc.tf.lower = S.bounds.t_final.lower;
bounds_nsc.tf.upper = S.bounds.t_final.upper;


%% Muscle-tendon forces
bounds_nsc.FTtilde.lower = zeros(1,NMuscle);
bounds_nsc.FTtilde.upper = 5*ones(1,NMuscle);

%% Time derivative of muscle activations
bounds_nsc.vA.lower = (-1/100*ones(1,NMuscle))./(ones(1,NMuscle)*model_info.muscle_info.tdeact);
bounds_nsc.vA.upper = (1/100*ones(1,NMuscle))./(ones(1,NMuscle)*model_info.muscle_info.tact);

%% Time derivative of muscle-tendon forces
bounds_nsc.dFTtilde.lower = -1*ones(1,NMuscle);
bounds_nsc.dFTtilde.upper = 1*ones(1,NMuscle);

%% Torque actuator activations
bounds_nsc.a_a.lower = -ones(1,model_info.ExtFunIO.jointi.nq.torqAct);
bounds_nsc.a_a.upper = ones(1,model_info.ExtFunIO.jointi.nq.torqAct);

%% Torque actuator excitations
bounds_nsc.e_a.lower = -ones(1,model_info.ExtFunIO.jointi.nq.torqAct);
bounds_nsc.e_a.upper = ones(1,model_info.ExtFunIO.jointi.nq.torqAct);





end % end of function