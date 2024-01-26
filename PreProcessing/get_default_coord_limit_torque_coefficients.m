function [K_pass,theta_pass] = get_default_coord_limit_torque_coefficients(coord_names)
% --------------------------------------------------------------------------
% get_default_coord_limit_torque_coefficients
%   Returns the default coefficients that determine the  coordinate limit 
%   torque. This torque is expressed in function of coordinate value q:
%       Tau_pass = K_pass(1)*exp(K_pass(2)*(q - theta_pass(2))) + ...
%                  K_pass(3)*exp(K_pass(4)*(q - theta_pass(1)))
% 
%   If none of the specified coefficients match the coordinate name, this
%   function will look for a close match (e.g. assign "hip_flex" bounds to
%   "hip_flexion_r" coordinate). 
%
%   Coefficient values are taken from: 
%       Anderson III, Frank Clayton. A dynamic optimization solution for a 
%       complete cycle of normal gait. The University of Texas at Austin, 1999.
%
% INPUT:
%   - coord_names -
%   * cell array of coordinate names
%
% OUTPUT:
%   - K_pass -
%   * Table with 4 K_pass values for each input coordinate. 
% 
%   - theta_pass -
%   * Table with 2 theta_pass values for each input coordinate.
%
% Original author: Lars D'Hondt
% Original date: 12/April/2022
%
% Last edit by:
% Last edit date:
% --------------------------------------------------------------------------

%% Predefined 
K.hip_flexion = [-2.44 5.05 1.51 -21.88];
theta.hip_flexion = [-0.6981 1.81];
K.hip_adduction = [-0.03 14.94 0.03 -14.94];
theta.hip_adduction = [-0.5 0.5];
K.hip_rotation = [-0.03 14.94 0.03 -14.94];
theta.hip_rotation = [-0.92 0.92];
K.knee_angle = [-6.09 33.94 11.03 -11.33];
theta.knee_angle = [-2.4 0.13];
K.ankle_angle = [-2.03 38.11 0.18 -12.12];
theta.ankle_angle = [-0.74 0.52];
K.subtalar_angle = [-60.21 16.32 60.21 -16.32];
theta.subtalar_angle = [-0.65 0.65];
K.mtp_angle = [-0.9 14.87 0.18 -70.08];
theta.mtp_angle = [0 65/180*pi];
K.lumbar_extension = [-0.35 30.72 0.25 -20.36];
theta.lumbar_extension = [-0.5235987755982988 0.17];
K.lumbar_bending = [-0.25 20.36 0.25 -20.36];
theta.lumbar_bending = [-0.3490658503988659 0.3490658503988659];
K.lumbar_rotation = [-0.25 20.36 0.25 -20.36];
theta.lumbar_rotation = [-0.3490658503988659 0.3490658503988659];

K.mtj_angle = [-50 15 60 -30]';
theta.mtj_angle = [-0.4 0.5]';



%% Match bounds to input coordinates
% names of coordinates with predefined bounds
coeff_names = fieldnames(K);

% default value is 0
K_pass = nan(length(coord_names),4);
theta_pass = nan(length(coord_names),2);

% loop through input coordinates
for i=1:length(coord_names)
    coord_name_i = coord_names{i};
    % In order of priority:
    % 1) exact match
    idx = find(strcmp(coeff_names,coord_names{i}));

    % 2) assume symmetry
    if isempty(idx)
        coord_name_i_lr = coord_name_i;
        if strcmp(coord_name_i_lr(end-1:end),'_l') || strcmp(coord_name_i_lr(end-1:end),'_r')
            coord_name_i_lr = coord_name_i_lr(1:end-2);
        end
        idx = find(contains(coeff_names(:),coord_name_i_lr));
    end

    % 3) bound is indicated by shortened notation
    if isempty(idx)
        for j=1:length(coeff_names)
            if contains(coord_name_i,coeff_names{j})
                idx = j;
                break
            end
        end
    end

    % Overwrite default values if specific bounds are found
    if ~isempty(idx)
        K_pass(i,:) = K.(coeff_names{idx(1)});
        theta_pass(i,:) = theta.(coeff_names{idx(1)});
    end
end




