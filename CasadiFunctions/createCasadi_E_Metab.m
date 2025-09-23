function [fgetMetabolicEnergySmooth2004all] = createCasadi_E_Metab(S,model_info)
% --------------------------------------------------------------------------
% createCasadi_E_Metab
%   Function to create Casadi functions for metabolic energy.
%   
% INPUT:
%   - S -
%   * setting structure S
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - fgetMetabolicEnergySmooth2004all -
%   * Casadi functions for metabolic energy based on
%   L. J. Bhargava, M. G. Pandy, and F. C. Anderson, "A phenomenological model
%   for estimating metabolic energy consumption in muscle contraction,” Journal 
%   of biomechanics, vol. 37, no. 1, pp. 81–88, 2004.
% 
% Original author: Ines Vandekerckhove, KU Leuven
% Original date: 01-12-2021
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.*
NMuscle = model_info.muscle_info.NMuscle;
FMo = struct_array_to_double_array(model_info.muscle_info.parameters,'FMo');
musclemass = struct_array_to_double_array(model_info.muscle_info.parameters,'muscle_mass');
pctst = struct_array_to_double_array(model_info.muscle_info.parameters,'slow_twitch_fiber_ratio');

%% Metabolic energy models
act_SX          = SX.sym('act_SX',NMuscle,1); % Muscle activations
exc_SX          = SX.sym('exc_SX',NMuscle,1); % Muscle excitations
lMtilde_SX      = SX.sym('lMtilde_SX',NMuscle,1); % N muscle fiber lengths
vM_SX           = SX.sym('vM_SX',NMuscle,1); % Muscle fiber velocities
Fce_SX          = SX.sym('FT_SX',NMuscle,1); % Contractile element forces
Fpass_SX        = SX.sym('FT_SX',NMuscle,1); % Passive element forces
Fiso_SX         = SX.sym('Fiso_SX',NMuscle,1); % N forces (F-L curve)

% Bhargava et al. (2004)
[energy_total_sm_SX,Adot_sm_SX,Mdot_sm_SX,Sdot_sm_SX,Wdot_sm_SX,...
    energy_model_sm_SX] = getMetabolicEnergySmooth2004all(exc_SX,act_SX,...
    lMtilde_SX,vM_SX,Fce_SX,Fpass_SX,musclemass,pctst,Fiso_SX,...
    FMo,model_info.mass,S.metabolicE.tanh_b,...
    [model_info.muscle_info.parameters.muscle_strength]');
fgetMetabolicEnergySmooth2004all = ...
    Function('fgetMetabolicEnergySmooth2004all',...
    {exc_SX,act_SX,lMtilde_SX,vM_SX,Fce_SX,Fpass_SX,Fiso_SX},...
    {energy_total_sm_SX,...
    Adot_sm_SX,Mdot_sm_SX,Sdot_sm_SX,Wdot_sm_SX,energy_model_sm_SX});

end