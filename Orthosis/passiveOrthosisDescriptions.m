function [fun_passiveOrthosis] = passiveOrthosisDescriptions(S_orthosis_settings_i)
% --------------------------------------------------------------------------
% passiveOrthosisDescriptions
%   This function contains functions describing different orthoses. When
%   called, this function selects the desired orthosis description and
%   returns its function handle.
%   Custom orthosis descriptions can be added as functions at the end of
%   this file. 
%   
% 
% INPUT:
%   - S_orthosis_settings_i -
%   * setting structure for the orthosis. "S_orthosis_settings_i.name" should
%   contain the name of the function you want to call when adding an
%   orthosis to your simulation model. 
%   To add an orthosis, in main.m, set "S.orthosis.settings{i}.name" to the 
%   name of your desired orthosis.
%   Add other fields than "name" to pass settings (e.g. stiffness
%   parameters) to your orthosis description.
%
%
% OUTPUT:
%   - fun_passiveOrthosis -
%   * handle of the function "settings_orthosis.name"
% 
% Original author: Lars D'Hondt
% Original date: 9/September/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

% get handle of the function corresponding to given string
fun_passiveOrthosis = str2func(S_orthosis_settings_i.name);


%% Add a function describing your orthosis below, or make a separate file in the same folder
% To be compatible, your orthosis description needs to be consistent with
% the following rules:
% INPUT:
%   - settings_orthosis -
%   * "S_orthosis__settings_i" will be passed
%
%   - q_ (coordinate name) -
%   * position of coordinate with this name. (m or rad)
%
%   - qdot_ (coordinate name) -
%   * velocity of coordinate with this name. ((m/s or rad/s)
%
% OUTPUT:
%   - T_ (coordinate name) -
%   * torque (or force) applied to this coordinate by the orthosis
%

% AFO with linear stiffness for right foot
function [T_ankle_angle_r,T_mtp_angle_r] = AFO_passive_r(settings_orthosis,q_ankle_angle_r,q_mtp_angle_r)
    k_ankle = settings_orthosis.ankle_stiffness; % ankle stiffness in Nm/rad
    k_mtp = settings_orthosis.mtp_stiffness; % mtp stiffness in Nm/rad

    T_ankle_angle_r = -k_ankle*q_ankle_angle_r;
    T_mtp_angle_r = -k_mtp*q_mtp_angle_r;
end

% AFO with linear stiffness for left foot
function [T_ankle_angle_l,T_mtp_angle_l] = AFO_passive_l(settings_orthosis,q_ankle_angle_l,q_mtp_angle_l)
    k_ankle = settings_orthosis.ankle_stiffness; % ankle stiffness in Nm/rad
    k_mtp = settings_orthosis.mtp_stiffness; % mtp stiffness in Nm/rad

    T_ankle_angle_l = -k_ankle*q_ankle_angle_l;
    T_mtp_angle_l = -k_mtp*q_mtp_angle_l;
end

% ...

end % end of passiveOrthosisDescriptions