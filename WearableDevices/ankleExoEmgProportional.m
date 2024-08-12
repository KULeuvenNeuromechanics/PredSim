function [exo] = ankleExoEmgProportional(init, settings_orthosis)
% --------------------------------------------------------------------------
% ankleExoEmgProportional
%   Ankle exoskeleton that applies a torque proportional to soleus
%   activation.
%
%
% INPUT:
%   - init -
%   * struct with information used to initialise the Orthosis object.
% 
%   - settings_orthosis -
%   * struct with information about this orthosis, containing the fields:
%       - function_name = 'ankleExoEmgProportional'  i.e. name of this function
%       - gain:             T = gain*act
%       - left_right:       'l' for left or 'r' for right
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

T_ankle = (emg-0.05)*gain; % only use activation above lower bound (0.05)
exo.addCoordForce(T_ankle,['ankle_angle_',side]);

end