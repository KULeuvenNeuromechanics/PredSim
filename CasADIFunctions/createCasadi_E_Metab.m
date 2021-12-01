function [fgetMetabolicEnergySmooth2004all] = createCasadi_E_Metab(NMuscle,MTparameters_m)
%% createCasadi_E_Metab.m
%Function to create Casadi functions for metabolic energy.
%
%INPUT -> need to be renamed
% - NMuscle
% - MTparameters_m
%
%OUTPUT
%  - Casadi function
%
% Authors: Ines Vandekerckhove, KU Leuven
% Date: 01-12-2021 

%% Metabolic energy models
act_SX          = SX.sym('act_SX',NMuscle,1); % Muscle activations
exc_SX          = SX.sym('exc_SX',NMuscle,1); % Muscle excitations
lMtilde_SX      = SX.sym('lMtilde_SX',NMuscle,1); % N muscle fiber lengths
vM_SX           = SX.sym('vM_SX',NMuscle,1); % Muscle fiber velocities
Fce_SX          = SX.sym('FT_SX',NMuscle,1); % Contractile element forces
Fpass_SX        = SX.sym('FT_SX',NMuscle,1); % Passive element forces
Fiso_SX         = SX.sym('Fiso_SX',NMuscle,1); % N forces (F-L curve)
musclemass_SX   = SX.sym('musclemass_SX',NMuscle,1); % Muscle mass
pctst_SX        = SX.sym('pctst_SX',NMuscle,1); % Slow twitch ratio
modelmass_SX    = SX.sym('modelmass_SX',1); % Model mass
b_SX            = SX.sym('b_SX',1); % Parameter determining tanh smoothness
% Bhargava et al. (2004)
[energy_total_sm_SX,Adot_sm_SX,Mdot_sm_SX,Sdot_sm_SX,Wdot_sm_SX,...
    energy_model_sm_SX] = getMetabolicEnergySmooth2004all(exc_SX,act_SX,...
    lMtilde_SX,vM_SX,Fce_SX,Fpass_SX,musclemass_SX,pctst_SX,Fiso_SX,...
    MTparameters_m(1,:)',modelmass_SX,b_SX);
fgetMetabolicEnergySmooth2004all = ...
    Function('fgetMetabolicEnergySmooth2004all',...
    {exc_SX,act_SX,lMtilde_SX,vM_SX,Fce_SX,Fpass_SX,musclemass_SX,...
    pctst_SX,Fiso_SX,modelmass_SX,b_SX},{energy_total_sm_SX,...
    Adot_sm_SX,Mdot_sm_SX,Sdot_sm_SX,Wdot_sm_SX,energy_model_sm_SX});

end