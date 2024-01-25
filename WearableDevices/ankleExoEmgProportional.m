function [exo] = ankleExoEmgProportional(init, settings_orthosis)
% --------------------------------------------------------------------------
% ankleExoZhang2017
%   Ankle exoskeleton that applies a torque profile (torque in function of
%   stride) to the ankle. 
%
%   This function requires additional dependencies, which can be downloaded
%   from: 
%   https://www.science.org/doi/10.1126/science.aal5054#supplementary-materials
%
%   References
%   [1] J. Zhang et al., “Human-in-the-loop optimization of exoskeleton 
%   assistance during walking,” Science, vol. 356, pp. 1280–1283, Jun. 2017, 
%   doi: 10.1126/science.aal5054.
%
% INPUT:
%   - init -
%   * struct with information used to initialise the Orthosis object.
% 
%   - settings_orthosis -
%   * struct with information about this orthosis, containing the fields:
%       - function_name = ankleExoZhang2017  i.e. name of this function   
%       - dependencies_path path to dependencies
%       - isFullGaitCycle   assistance profile for full stride when true,
%       half stride when false. Default is false.
%       - peak_torque:      peak torque in Nm/rad
%       - peak_timing:      timing of peak as % of stride
%       - rise_time:        rise time as % of stride
%       - drop_time:        drop time as % of stride
%   Values are set via S.orthosis.settings{i} in main.m
%
%
% OUTPUT:
%   - exo -
%   * an object of the class Orthosis
% 
% Original author: Lars D'Hondt
% Original date: 8/January/2024
% --------------------------------------------------------------------------

% create Orthosis object
exo = Orthosis('exo',init);


gain = settings_orthosis.gain;
side = settings_orthosis.left_right; % 'l' for left or 'r' for right

emg = exo.var_muscle(['soleus_',side]);

T_ankle = [0;0; (emg-0.05)*gain];

% apply exo torque on tibia and calcn
exo.addBodyMoment(T_ankle, ['T_exo_shank_',side],['tibia_',side]);
exo.addBodyMoment(-T_ankle, ['T_exo_foot_',side],['calcn_',side],['tibia_',side]);


end