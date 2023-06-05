function [Qs] = generate_dummy_motion(S,model_info,n_data_points)
% --------------------------------------------------------------------------
% generate_dummy_motion
%   Generate a set of coordinate values for an OpenSim model. Coordinate
%   values are taken in the range of motion. Upper and lower bounds are
%   determined from (in descending order of priority)
%       1) user input of individual bounds (S.misc.msk_geom_bounds)
%       2) user input of table with default bounds (S.misc.default_msk_geom_bounds)
%       3) coordinate value range from .osim file
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
%   - n_data_points - 
%   * number of data points (discrete time points)
%
% OUTPUT:
%   - Qs -
%   * Array with coordinate values with a column for each coordinate and a
%   row for each data point. Translations are in meters, rotations in
%   radians.
%
% Original author: Lars D'Hondt
% Original date: 05/April/2022
%
% Last edit by: Lars D'Hondt (Provide multiple ways to set bounds)
% Last edit date: 5/June/2023
% --------------------------------------------------------------------------

%% Set bounds for joint range of motion
% Initialise
Q_bounds = nan(2,model_info.ExtFunIO.jointi.nq.all);
Q_bounds(:,model_info.ExtFunIO.jointi.floating_base) = 0;

% Default upper and lower bounds of dummy motion
if exist(S.misc.default_msk_geom_bounds,'file')
    default_bounds = readtable(S.misc.default_msk_geom_bounds);

    for j=1:model_info.ExtFunIO.jointi.nq.all

        % search for coordinate name in the provided table with defaults
        default_bounds_j = default_bounds(strcmp(default_bounds.name,...
            model_info.ExtFunIO.coord_names.all{j}),:);
        if ~isempty(default_bounds_j)
            Q_bounds(1,j) = default_bounds_j.lower;
            Q_bounds(2,j) = default_bounds_j.upper;
        end

    end
end

bounds_def = Q_bounds';

% Any bounds that were not found in the table (i.e. values are NaN) are
% taken from the .osim file. 
import org.opensim.modeling.*;
model = Model(model_info.osim_path);

for j=1:model_info.ExtFunIO.jointi.nq.all
    coord_j = model.getCoordinateSet().get(model_info.ExtFunIO.coord_names.all{j});
    lb_j = coord_j.getRangeMin();
    ub_j = coord_j.getRangeMax();

    if any(model_info.ExtFunIO.jointi.rotations(:)==j)
        lb_j = lb_j*180/pi;
        ub_j = ub_j*180/pi;
    end
    if isnan(Q_bounds(1,j))
        Q_bounds(1,j) = lb_j;
    end
    if isnan(Q_bounds(2,j))
        Q_bounds(2,j) = ub_j;
    end
end

bounds_osim = Q_bounds';

% User inputs of individual bounds have highest priority, so we overwrite
% any previous value.
if ~isempty(S.misc.msk_geom_bounds)
    [new_lb,new_ub] = unpack_name_value_combinations(S.misc.msk_geom_bounds,...
        model_info.ExtFunIO.coord_names.all,[1,1]);
    for i=1:model_info.ExtFunIO.jointi.nq.all
        if ~isnan(new_lb(i))
            Q_bounds(1,i) = new_lb(i);
        end
        if ~isnan(new_ub(i))
            Q_bounds(2,i) = new_ub(i);
        end

    end
end

%% Create random poses
% Construct scale from bounds
Q_scale = diff(Q_bounds);

% Generate array with random coordinate values (range 0-1) for half the data points
n_data_points = round(n_data_points/2);
Qs1 = lhsdesign(n_data_points,model_info.ExtFunIO.jointi.nq.all);

% scale and offset random joint angles to fit bounds
Qs1 = Qs1.*(ones(n_data_points,1)*Q_scale) + ones(n_data_points,1)*Q_bounds(1,:);

% base does not have to move for analysis
Qs1(:,model_info.ExtFunIO.jointi.floating_base) = 0;

% angles from degrees to radians
Qs1(:,model_info.ExtFunIO.jointi.rotations) = Qs1(:,model_info.ExtFunIO.jointi.rotations)*pi/180;

% Create second array that is inverse/opposite of first array
Qs2 = zeros(size(Qs1));
Qs2(:,model_info.ExtFunIO.symQs.QsInvA) = Qs1(:,model_info.ExtFunIO.symQs.QsInvB);
Qs2(:,model_info.ExtFunIO.symQs.QsOpp) = -Qs1(:,model_info.ExtFunIO.symQs.QsOpp);

% total array of coordinate values contains symmetric entries
Qs = [Qs1; Qs2];



end % end of function