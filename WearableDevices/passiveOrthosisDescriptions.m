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

end % end of passiveOrthosisDescriptions

%% Add a function describing your orthosis below, or make a separate file in the same folder
% To be compatible, your orthosis description needs to be consistent with
% the following rules:
% INPUT:
%   - settings_orthosis -
%   * "S_orthosis__settings_i" will be passed
%
%   - q_(coordinate) -
%   * position of coordinate with this name. (m or rad)
%
%   - qdot_(coordinate) -
%   * velocity of coordinate with this name. (m/s or rad/s)
%
% OUTPUT:
%   - input_names -
%   * cell array containing strings q_(coordinate) and qdot_(coordinate) in the
%   order they appear as inputs. Entries should match coordinate names from
%   the selected OpenSim model exactly.
%
%   - output_names -
%   * cell array containing strings (coordinate) in the order they appear as 
%   outputs. Entries should match coordinate names from the selected OpenSim
%   model exactly.
%
%   - T_(coordinate) -
%   * torque (or force) applied to this coordinate by the orthosis
%

% AFO with linear stiffness
function [input_names,output_names,T_ankle_angle,T_mtp_angle] = AFO_passive(settings_orthosis,q_ankle_angle,q_mtp_angle)

% read settings
    % ankle stiffness in Nm/rad
    k_ankle = settings_orthosis.ankle_stiffness;
    % mtp stiffness in Nm/rad
    k_mtp = settings_orthosis.mtp_stiffness; 
    % left or right
    side = settings_orthosis.left_right;

% calculate outputs
    T_ankle_angle = -k_ankle*q_ankle_angle;
    T_mtp_angle = -k_mtp*q_mtp_angle;

% input and output names
    input_names = {['q_ankle_angle_' side],['q_mtp_angle_' side]};
    output_names = {['ankle_angle_' side],['mtp_angle_' side]};

end % end of AFO_passive

% AFO with linear stiffness, toggled via body weight clutch
function [input_names,output_names,T_ankle_angle,T_mtp_angle] = AFO_passive_BWC(settings_orthosis,q_ankle_angle,q_mtp_angle,GRF_y)

% read settings
    % lower threshold to activate clutch
    lower_threshold = settings_orthosis.lower_threshold;
    % upper threshold to activate clutch
    upper_threshold = settings_orthosis.upper_threshold;
    % left or right
    side = settings_orthosis.left_right;

% calculate outputs
    % toggle clutch
    tglcl = smoothIf(GRF_y,lower_threshold,upper_threshold);
    % spring force without clutch
    [input_names,output_names,T_ankle_angle,T_mtp_angle] = AFO_passive(settings_orthosis,q_ankle_angle,q_mtp_angle);
    % apply clutch
    T_ankle_angle = T_ankle_angle*tglcl;

% input and output names
    input_names{end+1} = ['GRF_y_' side];

end % end of AFO_passive_BWC


% add functions here ...

%% Helper functions
function y = smoothIf(x,lower_threshold,upper_threshold)
% Reformulate conditional statement to be compatible with algorithmic
% differentiation.
%
% if x <= lower_threshold
%   y = 0; % false
% elseif x >= upper_threshold
%   y = 1; % true
% else
%   y > 0 & y < 1; % transient
% end
%
% take lower_threshold > upper_threshold to flip the true and false outputs
%
%--------------------------------------------------------------------------
    % range that is affected by smoothing
    range = upper_threshold - lower_threshold;

    % midpoint of  range
    mid = (lower_threshold + upper_threshold)/2;

    % scale x to range
    x_scaled = (x - mid)/range + 1/2;

    % apply tanh-smoothing
    y = (tanh((2*x_scaled-1)*pi)+1)/2;

end
