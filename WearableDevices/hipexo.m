function [exo] = hipexo(init, settings_orthosis)
% --------------------------------------------------------------------------
% Provisional hip exo with 2 hip flexion/extension actuators
% Example use of adding internal dynamics to the orthosis 
% using addDynamics method, and adding controls to optimize for (var_opti
% method)
%
% INPUT:
%   - init -
%   * struct with information used to initialise the Orthosis object.
% 
%   - settings_orthosis -
%   * struct with information about this orthosis, containing the fields:
%       - function_name = hipexo  i.e. name of this function   
%       - gain: EMG gain
%       - left_right: char, 'l' or 'r'
%       - dynamics: struct with lower and upper limits for opti controls and
%           states
%   Values are set via S.orthosis.settings{i} in main.m, with i the index
%   of the orthosis.
%
%
% OUTPUT:
%   - exo -
%   * an object of the class Orthosis
% 
% Original author: Sander De Groof
% Original date: 3/October/2025
% --------------------------------------------------------------------------

% create Orthosis object
exo = Orthosis('exo',init, true);

% read settings that were passed from main.m
if isfield(settings_orthosis,'isFullGaitCycle')
    isFullGaitCycle = settings_orthosis.isFullGaitCycle;
else
    isFullGaitCycle = false;
end

side = settings_orthosis.left_right; % 'l' for left or 'r' for right
% gain = settings_orthosis.gain;
% 
% emg1 = exo.var_muscle(['glut_max1_',side]);
% emg2 = exo.var_muscle(['glut_max2_',side]);
% emg3 = exo.var_muscle(['glut_max3_',side]);
% 
% emg = (emg1 + emg2 + emg3)/3;

%u_hipfl_emg = -(emg-0.05)*gain; % only use activation above lower bound (0.05)

%test time variable control
% timebase = 0:init.Nmesh;
% if strcmp(side,'l')
%     TimeVarControl = 10*sin(2*pi*timebase/timebase(end)); %10Nm torque in the middle of cycle
% else
%     TimeVarControl = 10*sin(2*pi*timebase/timebase(end)+pi); %10Nm torque in the middle of cycle
% end

% Simple first order dynamics dx/dt = (u-x)/tau
tc_tau = 0.05; % time constant in seconds
state_x = exo.var_opti(['state_x_' side '_side'],'state',[settings_orthosis.dynamics.xl settings_orthosis.dynamics.xu]);
state_x2 = exo.var_opti(['state_x2_' side '_side'],'state',[settings_orthosis.dynamics.xl settings_orthosis.dynamics.xu]);
control_u = exo.var_opti(['control_u_' side '_side'],'control',[settings_orthosis.dynamics.ul settings_orthosis.dynamics.uu]);

%control_exo = u_hipfl_emg + control_u; % to make it interesting
control_exo = control_u;

%exo.addDynamics((control_exo-state_x)/tc_tau,['state_x_' side '_side']); %state and control
exo.addDynamics([(control_exo-state_x)/tc_tau;(control_exo-state_x2)/tc_tau],{['state_x_' side '_side'],['state_x2_' side '_side']}); %state and control
%exo.addDynamics((emg-state_x)/tc_tau,['state_x_' side '_side']); %state, but no control
exo.addCoordForce(state_x,['hip_flexion_',side]);
%exo.addCoordForce(state_x+TimeVarControl,['hip_flexion_',side]);
%exo.addCoordForce(control_u,['hip_flexion_',side]); % add control straight away

%exo.addVarToPostProcessing(state_x,['state_x_' side '_side'])
%exo.addVarToPostProcessing(control_exo,['control_u_' side '_side'])
%exo.addVarToPostProcessing(control_u,['control_u_' side '_side'])
end

