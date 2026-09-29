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
%   Two sets of coordinate headers are stored:
%
%       R.colheaders.coordinates
%           -> headers for the KINEMATIC curves (R.kinematics.Qs)
%
%       R.colheaders.coordinates_kinet
%           -> headers for the KINETIC curves (R.kinetics.T_ID)
%
%   Both header sets initially contain the original PredSim coordinate
%   names. After the transformations, the names are changed according
%   to the corresponding clinical convention.
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
% Revised date: September 02, 2026

%%

    %% Default = right side
    if nargin < 2 || isempty(side)
        side = 'right';
    end

    % Check side
    if ~strcmpi(side, 'right') && ~strcmpi(side, 'left')
        error('side must be either ''right'' or ''left''.');
    end


    %% Load data
    data = load(result_path);
    R = data.R;


    % Helper: find columns whose coordinate name contains any of
    % the given substrings
    findCols = @(substrings) any(cell2mat(cellfun( ...
        @(s) contains(R.colheaders.coordinates, s), ...
        substrings, ...
        'UniformOutput', false)), 2);

    %%  KINEMATICS
    %% Flip signs of selected kinematics

    % Always flip these
    coords_kin_flip = { ...
        'pelvis_tilt', ...
        'knee', ...
        'lumbar_extension'};

    % Side-dependent coordinates
    if strcmpi(side, 'right')

        coords_kin_flip = [coords_kin_flip, { ...
            'pelvis_list', ...
            'lumbar_bending'}];

    elseif strcmpi(side, 'left')

        coords_kin_flip = [coords_kin_flip, { ...
            'pelvis_rotation', ...
            'lumbar_rotation'}];

    end

    idx = findCols(coords_kin_flip);

    R.kinematics.Qs(:,idx) = ...
        -R.kinematics.Qs(:,idx);


    %% Apply 13 degree offsets

    % Add 13° to pelvis tilt and hip flexion
    coords = {'pelvis_tilt', 'hip_flexion'};

    idx = findCols(coords);

    R.kinematics.Qs(:,idx) = ...
        R.kinematics.Qs(:,idx) + 13;


    % Subtract 13° from lumbar extension
    coords = {'lumbar_extension'};

    idx = findCols(coords);

    R.kinematics.Qs(:,idx) = ...
        R.kinematics.Qs(:,idx) - 13;


    %%  KINETICS
    %% Flip signs of selected kinetics

    coords_kinet_flip = { ...
        'hip_flexion', ...
        'hip_adduction', ...
        'ankle'};

    idx = findCols(coords_kinet_flip);

    R.kinetics.T_ID(:,idx) = ...
        -R.kinetics.T_ID(:,idx);


    % Make a copy for kinetics before changing any kinematic names.

    R.colheaders.coordinates_kinet = ...
        R.colheaders.coordinates;


    %%  RENAME KINEMATIC COORDINATES

    rename_kinematics = {
        'pelvis_tilt',       'pelvis_ant_tilt';
        'knee_angle',        'knee_flexion';
        'lumbar_extension',  'lumbar_flexion';
        'pelvis_rotation',   'pelvis_int_rotation';
        'lumbar_rotation',   'lumbar_int_rotation';
        'pelvis_list',       'pelvis_up';
        'lumbar_bending',    'lumbar_up'
    };


    for i = 1:size(rename_kinematics, 1)

        old_name = rename_kinematics{i,1};
        new_name = rename_kinematics{i,2};

        % Only rename coordinates that were actually sign-flipped
        % for the selected side
        if any(cellfun(@(x) contains(x, old_name) || ...
                            contains(old_name, x), ...
                            coords_kin_flip))

            R.colheaders.coordinates = cellfun( ...
                @(x) strrep(x, old_name, new_name), ...
                R.colheaders.coordinates, ...
                'UniformOutput', false);

        end
    end


    %%  RENAME KINETIC COORDINATES

    rename_kinetics = {
        'hip_flexion',    'hip_extension';
        'hip_adduction',  'hip_abduction';
        'ankle',          'ankle_plantar_flexion'
    };


    for i = 1:size(rename_kinetics, 1)

        old_name = rename_kinetics{i,1};
        new_name = rename_kinetics{i,2};

        R.colheaders.coordinates_kinet = cellfun( ...
            @(x) strrep(x, old_name, new_name), ...
            R.colheaders.coordinates_kinet, ...
            'UniformOutput', false);

    end


    %% Put transformed R back into data

    data.R = R;

    %%  SAVE

    [path, name, ext] = fileparts(result_path);

    if strcmpi(side, 'right')
        suffix = '_transf_right';
    else
        suffix = '_transf_left';
    end

    new_path = fullfile(path, [name suffix ext]);

    save(new_path, '-struct', 'data');

end