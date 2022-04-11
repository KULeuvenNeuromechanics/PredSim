function [Qs] = generate_dummy_motion(S,model_info,n_data_points)
% --------------------------------------------------------------------------
% generate_dummy_motion
%   Generate a set of coordinate values for an OpenSim model.
% 
% INPUT:
%   - S -
%   * setting structure S
%
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
%   - n_data_points - 
%   * number of data points (discrete time points)
%
% OUTPUT:
%   - Qs -
%   * Table with coordinate values with a column for each coordinate and a
%   row for each data point. Translations are in meters, rotations in
%   radians.
% 
% note: The boundaries of the dummy motion are hard-coded. They should
%       become user inputs in a future update;
%
% Original author: Lars D'Hondt
% Original date: 05/April/2022
%
% Last edit by:
% Last edit date:
% --------------------------------------------------------------------------



%% Create dummy motion

% bounds of training dataset (in degrees)
% hard-coded for now, make input later
Bounds(1,:) = {'hip_flex',[-50 50]};
Bounds(2,:) = {'hip_add',[-30 30]};
Bounds(3,:) = {'hip_rot',[-30 30]};
Bounds(4,:) = {'knee',[-90 0]};
Bounds(5,:) = {'ankle',[-30 30]};
Bounds(6,:) = {'subt',[-30 30]};
Bounds(7,:) = {'mtj',[-30 30]};
Bounds(8,:) = {'mtp',[-20 50]};
Bounds(9,:) = {'lumbar_ext',[-30 30]};
Bounds(10,:) = {'lumbar_bend',[-30 30]};
Bounds(11,:) = {'lumbar_rot',[-30 30]};


% names
coordinate_names = model_info.ExtFunIO.coord_names.all;

% sizes and indices
n_coord = length(coordinate_names);

% default bounds
Q_bounds = [-30;30]*ones(1,n_coord);

% adjust bounds
for j=1:size(Bounds,1)
    idx = find(contains(coordinate_names,Bounds{j,1}));
    if ~isempty(idx)
        Q_bounds(:,idx) = reshape(Bounds{j,2},2,1)*ones(1,length(idx));
    end
end

% construct scale from bounds
Q_scale = diff(Q_bounds);

% generate random joint angles (range 0-1)
Qs = lhsdesign(n_data_points,n_coord);


% scale and offset random joint angles to fit bounds
Qs = Qs.*(ones(n_data_points,1)*Q_scale) + ones(n_data_points,1)*Q_bounds(1,:);


% base does not have to move for analysis
Qs(:,model_info.ExtFunIO.jointi.floating_base) = 0;

% angles from degrees to radians
Qs(:,model_info.ExtFunIO.jointi.rotations) = Qs(:,model_info.ExtFunIO.jointi.rotations)*pi/180;


