function [coord_bounds] = get_bounds_dummy_motion(coord_names)
% --------------------------------------------------------------------------
% get_bounds_dummy_motion
%   Returns the upper and lower bound of coordinates. These bounds are used 
%   to generate a dummy motion to determine the musculoskeletal geometry.
%
%   THESE ARE NOT THE BOUNDS USED TO FORMULATE AND SCALE THE OCP!!!!
% 
%   If none of the specified bounds match the coordinate name, this
%   function will look for a close match (e.g. assign "hip_flex" bounds to
%   "hip_flexion_r" coordinate). If no close enough match is found, the
%   bounds are -30 and 30;
%
% INPUT:
%   - coord_names -
%   * cell array of coordinate names
%
% OUTPUT:
%   - bounds -
%   * Table with lower and upper bound for each input coordinate. Translations
%   are in meters, rotations in radians.
% 
%
% Original author: Lars D'Hondt
% Original date: 14/April/2022
%
% Last edit by:
% Last edit date:
% --------------------------------------------------------------------------

%% Predefined bounds
% Template:
% Bounds.(coordinate name) = [lower bound, upper bound];
Bounds.hip_flex = [-50 50];
Bounds.hip_add = [-30 30];
Bounds.hip_rot = [-30 30];
Bounds.knee = [-90 0];
Bounds.ankle = [-30 30];
Bounds.subt = [-30 30];
Bounds.mtj = [-30 30];
Bounds.mtp = [-20 50];
Bounds.lumbar_ext = [-30 30];
Bounds.lumbar_bend = [-30 30];
Bounds.lumbar_rot = [-30 30];
% ... add more fields if needed



%% Match bounds to input coordinates
% names of coordinates with predefined bounds
bound_names = fieldnames(Bounds);

% default value is 30
coord_bounds = [-30;30].*ones(1,length(coord_names));

% loop through input coordinates
for i=1:length(coord_names)
    coord_name_i = coord_names{i};
    % In order of priority:
    % 1) exact match
    idx = find(strcmp(bound_names,coord_names{i}));

    % 2) assume symmetry
    if isempty(idx)
        coord_name_i_lr = coord_name_i;
        if strcmp(coord_name_i_lr(end-1:end),'_l') || strcmp(coord_name_i_lr(end-1:end),'_r')
            coord_name_i_lr = coord_name_i_lr(1:end-2);
        end
        idx = find(contains(bound_names(:),coord_name_i_lr));
    end

    % 3) bound is indicated by shortened notation
    if isempty(idx)
        for j=1:length(bound_names)
            if contains(coord_name_i,bound_names{j})
                idx = j;
                break
            end
        end
    end

    % Overwrite default values if specific bounds are found
    if ~isempty(idx)
        coord_bounds(:,i) = Bounds.(bound_names{idx(1)});
    end
end




