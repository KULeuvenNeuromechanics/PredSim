function R = transform_result_clinical_format(result_path, side)
%   Transforms the kinematic and kinetic curves of a PredSim result R
%   (stored in a matfile) into the clinical format used by Vicon/Nexus.
%   The transformed result is saved to a new matfile with suffix
%   '_transf', alongside the unchanged original variables
%   (model_info, setup, stats, w_opt).
%
%   Sign flips:
%       - Kinematics (R.kinematics.Qs): pelvis_tilt, knee, and
%         lumbar_extension are always sign-flipped. Additionally,
%         pelvis_list and lumbar_bending are flipped if side is 'right',
%         or pelvis_rotation and lumbar_rotation are flipped if side is
%         'left'.
%       - Kinetics (R.kinetics.T_ID): hip_flexion, hip_adduction, and
%         ankle moments are sign-flipped, regardless of side.
%       - Matching is done with `contains`, so both left and right
%         variants of a coordinate (e.g. hip_flexion_l/hip_flexion_r)
%         are affected whenever the base name is listed above.
%
%   Offsets:
%       - 13° is added to pelvis_tilt and hip_flexion.
%       - 13° is subtracted from lumbar_extension.
%
% INPUT:
%   result_path - path to matfile with PredSim results
%   side        - 'right' (default, used if omitted or empty) or 'left'.
%                 Determines which side-specific coordinates get flipped.
%
% OUTPUT:
%   R - PredSim results transformed to clinical format (also returned
%       directly, in addition to being saved to disk)
%   A new matfile is saved at result_path with suffix '_transf', containing
%   the transformed R plus the original model_info, setup, stats, and w_opt
%   unchanged.
%
% Original author: Ines Vandekerckhove
% Original date: August 18, 2026

%%
    % Default = right side
    if nargin < 2 || isempty(side)
        side = 'right';
    end

    data = load(result_path);
    R = data.R;

    % Helper: find columns whose coordinate name contains any of the given substrings
    findCols = @(substrings) any(cell2mat(cellfun(@(s) contains(R.colheaders.coordinates, s), substrings, 'UniformOutput', false)), 2);

    %% Flip signs of selected kinematics

    % Always flip these
    coords = {'pelvis_tilt', 'knee', 'lumbar_extension'};

    % Side-dependent coordinates
    if strcmpi(side, 'right')
        coords = [coords, {'pelvis_list', 'lumbar_bending'}];

    elseif strcmpi(side, 'left')
        coords = [coords, {'pelvis_rotation', 'lumbar_rotation'}];

    else
        error('side must be either ''right'' or ''left''.');
    end

    idx = findCols(coords);
    R.kinematics.Qs(:,idx) = -R.kinematics.Qs(:,idx);


    %% Flip signs of selected kinetics

    coords = {'hip_flexion', 'hip_adduction', 'ankle'};
    idx = findCols(coords);
    R.kinetics.T_ID(:,idx) = -R.kinetics.T_ID(:,idx);


    %% Apply 13 degree offsets

    % Add 13° to pelvis tilt and hip flexion
    coords = {'pelvis_tilt', 'hip_flexion'};
    idx = findCols(coords);
    R.kinematics.Qs(:,idx) = R.kinematics.Qs(:,idx) + 13;

    % Subtract 13° from lumbar extension
    coords = {'lumbar_extension'};
    idx = findCols(coords);
    R.kinematics.Qs(:,idx) = R.kinematics.Qs(:,idx) - 13;

    data.R = R;

   %% Save to new file with '_transf' suffix
    [path, name, ext] = fileparts(result_path);

    if strcmpi(side, 'right')
        suffix = '_transf_right';
    elseif strcmpi(side, 'left')
        suffix = '_transf_left';
    else
        error('side must be either ''right'' or ''left''.');
    end

    new_path = fullfile(path, [name suffix ext]);
    save(new_path, '-struct', 'data');
    
end

 