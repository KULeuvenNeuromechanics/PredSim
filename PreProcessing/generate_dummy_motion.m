function [Qs] = generate_dummy_motion(S,model_info,n_data_points)
% --------------------------------------------------------------------------
% generate_dummy_motion
%   Generate a set of coordinate values for an OpenSim model.
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
% Last edit by:
% Last edit date:
% --------------------------------------------------------------------------


% Default upper and lower bounds of dummy motion
Q_bounds = get_default_bounds_dummy_motion(model_info.ExtFunIO.coord_names.all);

% adapt bounds based on user input
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



