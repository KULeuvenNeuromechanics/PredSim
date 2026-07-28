% This function computes the muscle energy expenditure based on the model
% of Margaria. 
%
% Author: Tim van der Zee
% Date: 28/07/2026
%
% INPUTS:
%   exc: muscle excitations
%   act: muscle activations
%   lMtilde: normalized muscle fiber lengths
%   vMtilde: normalized muscle fiber velocities (positive is lengthening)
%       Note: defined as vM/lMopt, which is different than typical vM/vMax
%   vM: muscle fiber velocities (positive is lengthening)
%   Fce: muscle force from contractile element 
%       (length + velocity components but not passive component; FMo.*Fcetilde)
%   musclemass: mass of muscle, mass = PCSA*rho*lMopt, PCSA = Fmax/sigma
%       Note: rho is density (1058.7 kg/m³) and sigma is tension
%   pctst: percentage of slow twitch fibers (0-1)
%   vcemax: maximal muscle fiber contraction velocities (default is 10)
%   Fiso: normalized muscle forces from active f-l relationship (FMltilde)
%   modelmass: mass of the musculoskeletal model
%   b: parameter determining transition smoothness for tanh approximations
%
% OUTPUTS:
%   energy_total: total metabolic energy rate
%   energy_am: energy rate from activation and maintenance
%   energy_sl: energy rate from shortening and lengthening
%   energy_mech: energy rate from mechanical work
%   energy_model: energy rate from energy_total including basal rate
%       energy_total: energy rate from totalHeatRate and energy_mech
%       totalHeatRate: energy rate from energy_am and energy_sl
%       Note: energy_total might be different than the sum of the different
%             components if the total heat rate was clamped to one.

function [energy_total,energy_a, energy_am,energy_sl,energy_mech,energy_model] = ...
    getMetabolicEnergyMargaria(exc,act,lMtilde,vM,lmo,Fce,Fpe, ...
        musclemass,pctst,vcemax,Fiso,modelmass,b,strength)
            
    
    % following getMetabolicEnergySmooth2004all
    musclemass=musclemass.*strength;
    
    %% Mechanical work
    % Include negative mechanical work
    energy_mech = -1*Fce.*vM./musclemass; % same as Umberger et al. (2003)

    energy_am = zeros(size(energy_mech));
    energy_sl = zeros(size(energy_mech));
    
    % based on 25% and -120% efficiency
    totalHeatRate = energy_mech(energy_mech > 0) * 3 + energy_mech(energy_mech < 0) * -1;
    
    %% Account for muscle mass
    energy_am = energy_am.*musclemass;
    energy_sl = energy_sl.*musclemass;
    energy_mech = energy_mech.*musclemass;
    energy_total = totalHeatRate.*musclemass+energy_mech;
    energy_a = zeros(size(energy_am)); % lumped into am

    %% Energy model
    basal_coef = 1.2; % default in OpenSim
    basal_exp = 1; % default in OpenSim
    energy_model = basal_coef*modelmass^basal_exp + sum(energy_total);

end
