function [exo] = ankleExoNuckols2020(init, settings_orthosis)
% --------------------------------------------------------------------------
% ankleExoNuckols2020
%   Defines a passive ankle exoskeleton as described by Nuckols et al.
% 
%   [1] R. W. Nuckols and G. S. Sawicki, “Impact of elastic ankle exoskeleton
%       stiffness on neuromechanics and energetics of human walking across 
%       multiple speeds,” Journal of NeuroEngineering and Rehabilitation, 
%       vol. 17, no. 1, p. 75, Jun. 2020, doi: 10.1186/s12984-020-00703-4.
%
%
% INPUT:
%   - init -
%   * struct with information used to initialise the Orthosis object.
% 
%   - settings_orthosis -
%   * struct with information about this orthosis, containing the fields:
%       - function_name = 'ankleExoNuckols2020'  i.e. name of this function           
%       - ankle_stiffness:  ankle stiffness in Nm/rad
%       - left_right:       'l' for left or 'r' for right
%   Values are set via S.orthosis.settings{i} in main.m
%
%
% OUTPUT:
%   - AFO -
%   * an object of the class Orthosis
% 
% Original author: Lars D'Hondt
% Original date: 5/January/2024
% --------------------------------------------------------------------------

% create Orthosis object
exo = Orthosis('exo',init);

% read settings that were passed from main.m
k_ankle = settings_orthosis.ankle_stiffness; % ankle stiffness in Nm/rad
side = settings_orthosis.left_right; % 'l' for left or 'r' for right

% get joint angles
q_ankle = exo.var_coord(['ankle_angle_',side]); % ankle angle in rad;

% calculate moments
T_ankle = k_ankle*(q_ankle+0.1).*smoothIf(q_ankle+0.1,0.05,0);
T_ankle = [0;0; T_ankle];

% apply exo torque on tibia and calcn
exo.addBodyMoment(T_ankle, ['T_exo_shank_',side],['tibia_',side]);
exo.addBodyMoment(-T_ankle, ['T_exo_foot_',side],['calcn_',side],['tibia_',side]);


end