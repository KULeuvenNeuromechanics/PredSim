function [exo] = EXO001_ankleExoNuckols2020(init, settings_orthosis)
% --------------------------------------------------------------------------
% parametricAFO
%   Defines an ankle-foot orthosis as a rotational spring on the ankle angle
%   and a rotational spring on the MTP angle. Each spring is defined by a
%   stiffness constant.
% 
%
% INPUT:
%   - init -
%   * struct with information used to initialise the Orthosis object.
% 
%   - settings_orthosis -
%   * struct with information about this orthosis, containing the fields:
%       - function_name = parametricAFO  i.e. name of this function           
%       - ankle_stiffness:  ankle stiffness in Nm/rad
%       - mtp_stiffness:    mtp stiffness in Nm/rad
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
exo = EXO001('exo',43,init);

% read settings that were passed from main.m
k_ankle = settings_orthosis.ankle_stiffness; % ankle stiffness in Nm/rad
side = settings_orthosis.left_right; % 'l' for left or 'r' for right
pf_offset = settings_orthosis.pf_offset; % neutral angle of spring in °
smooth_b = settings_orthosis.smooth_b; 
% exo.setFootForceBody('calcn_r');

% set side
exo.setSide(side);

% get joint angles
q_enc = exo.var_encoder("pos"); % plantarflexion is positive

% calculate moments
Deltaq = q_enc - pf_offset*pi/180;
Deltaq_pos = getSmoothingHuberCon(Deltaq,Deltaq,0,smooth_b,-1);
T_ankle = -k_ankle*Deltaq_pos;

% set moment
exo.setPlantarFlexionMoment(T_ankle);

end