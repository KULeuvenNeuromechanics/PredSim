function [bounds] = getBounds(S,model_info)
% --------------------------------------------------------------------------
% getBounds
%   This script provides bounds for the (unscaled) optimisation variables.
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - bounds -
%   * boundaries for all optimisation variables
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
    bounds.(i).lower = nan(1,NCoord);
    bounds.(i).upper = nan(1,NCoord);
end

%% Get default bounds
if exist(S.bounds.default_coordinate_bounds,'file')
    default_bounds = readtable(S.bounds.default_coordinate_bounds);
    for j=1:NCoord
        default_bounds_j = default_bounds(strcmp(default_bounds.name,coordinate_names{j}),:);
        if ~isempty(default_bounds_j)
            for i=coords
                bounds.(i).lower(j) = default_bounds_j.([char(i) '_lower']);
                bounds.(i).upper(j) = default_bounds_j.([char(i) '_upper']);    
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
    range_j = ub_j - lb_j;

    if isnan(bounds.Qs.lower(j))
        bounds.Qs.lower(j) = lb_j;
    end
    if isnan(bounds.Qs.upper(j))
        bounds.Qs.upper(j) = ub_j;
    end
    if isnan(bounds.Qdots.lower(j))
        bounds.Qdots.lower(j) = -range_j*10;
    end
    if isnan(bounds.Qdots.upper(j))
        bounds.Qdots.upper(j) = range_j*10;
    end
    if isnan(bounds.Qdotdots.lower(j))
        bounds.Qdotdots.lower(j) = -range_j*155;
    end
    if isnan(bounds.Qdotdots.upper(j))
        bounds.Qdotdots.upper(j) = range_j*155;
    end


end


%% Manual adjustments
% Pelvis_ty
bounds.Qs.upper(model_info.ExtFunIO.jointi.floating_base(5)) = model_info.IG_pelvis_y*1.2;
bounds.Qs.lower(model_info.ExtFunIO.jointi.floating_base(5)) = model_info.IG_pelvis_y*0.5;


% We adjust some bounds when we increase the speed to allow for the
% generation of running motions.
if S.subject.v_pelvis_x_trgt > 1.33
    % Shoulder flexion
    bounds.Qs.lower(idx_shoulder_flex) = -50*pi/180;
    % Pelvis tx
    bounds.Qdots.upper(model_info.ExtFunIO.jointi.base_forward) = ...
        bounds.Qdots.upper(model_info.ExtFunIO.jointi.base_forward)*2;
end

if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    bounds.Qs.upper(model_info.ExtFunIO.jointi.base_forward) = ...
        bounds.Qs.upper(model_info.ExtFunIO.jointi.base_forward)/2;
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
                bounds.(i).lower(coord_idx) = new_lb(j);
            end
        
            if ~isnan(new_ub(j))
                bounds.(i).upper(coord_idx) = new_ub(j);
            end
        end
    
    end
end



%% Hard bounds
% We impose the initial position of pelvis_tx to be 0
bounds.Qs_0.lower = bounds.Qs.lower;
bounds.Qs_0.upper = bounds.Qs.upper;
bounds.Qs_0.lower(model_info.ExtFunIO.jointi.base_forward) = 0;
bounds.Qs_0.upper(model_info.ExtFunIO.jointi.base_forward) = 0;
bounds.Qs_0.lower(model_info.ExtFunIO.jointi.base_lateral) = 0;
bounds.Qs_0.upper(model_info.ExtFunIO.jointi.base_lateral) = 0;



%% Muscle activations
bounds.a.lower = S.bounds.activation_all_muscles.lower*ones(1,NMuscle);
bounds.a.upper = ones(1,NMuscle);

if ~isempty(S.bounds.activation_selected_muscles)
    [new_lb,new_ub] = unpack_name_value_combinations(S.bounds.activation_selected_muscles,...
        muscle_info.muscle_names,[1,1]);
    for j=1:NMuscle
        if ~isnan(new_lb(j))
            bounds.a.lower(j) = new_lb(j);
        end
        if ~isnan(new_ub(j))
            bounds.a.upper(j) = new_ub(j);
        end
    end 
end

%% Muscle-tendon forces
bounds.FTtilde.lower = zeros(1,NMuscle);
bounds.FTtilde.upper = 5*ones(1,NMuscle);

%% Time derivative of muscle activations
bounds.vA.lower = (-1/100*ones(1,NMuscle))./(ones(1,NMuscle)*model_info.muscle_info.tdeact);
bounds.vA.upper = (1/100*ones(1,NMuscle))./(ones(1,NMuscle)*model_info.muscle_info.tact);

%% Time derivative of muscle-tendon forces
bounds.dFTtilde.lower = -1*ones(1,NMuscle);
bounds.dFTtilde.upper = 1*ones(1,NMuscle);

%% Torque actuator activations
bounds.a_a.lower = -ones(1,model_info.ExtFunIO.jointi.nq.torqAct);
bounds.a_a.upper = ones(1,model_info.ExtFunIO.jointi.nq.torqAct);

%% Torque actuator excitations
bounds.e_a.lower = -ones(1,model_info.ExtFunIO.jointi.nq.torqAct);
bounds.e_a.upper = ones(1,model_info.ExtFunIO.jointi.nq.torqAct);

%% Final time
bounds.tf.lower = S.bounds.t_final.lower;
if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    bounds.tf.upper = S.bounds.t_final.upper/2;
else
    bounds.tf.upper = S.bounds.t_final.upper;
end


end % end of function