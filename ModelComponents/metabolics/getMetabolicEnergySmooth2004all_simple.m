function [energy_total,Adot,Mdot,Sdot,Wdot,energy_model] = ...
    getMetabolicEnergySmooth2004all_simple(exc,act,lMtilde,vM,Fce,Fpass,...
        musclemass,pctst,Fiso,Fmax,modelmass,b,strength, lmo, vcemax)
% --------------------------------------------------------------------------
% getMetabolicEnergySmooth2004all_simple
%    This function computes the muscle energy expenditure based on the model
%    of Bhargava et al. (2004). This is a simpler form of the model that
%    only requires 1 tanh function.
%
%    Note: we did not include the length dependence of the maintenance heat rate.
%      In OpenSim, they use a Piecewise Linear Function
%      https://github.com/opensim-org/opensim-core/blob/master/OpenSim/Simulation/Model/Bhargava2004MuscleMetabolicsProbe.cpp#L103
%    Note: we did consider lengthening heat rate, i.e. 
%      get_use_force_dependent_shortening_prop_constant is true (OpenSim)
%   
% INPUT:
%   - exc -
%   * muscle excitations
%
%   - act -
%   * muscle activations
%
%   - lMtilde -
%   * normalized muscle fiber lengths
%
%   - vM -
%   * muscle fiber velocities (positive is lengthening)
%
%   - Fce -
%   * muscle force from contractile element 
%       (length + velocity components but not passive component; FMo.*Fcetilde)
%
%   - Fpass -
%   * muscle force from passive element (FMo.*Fpetilde)
%
%   - musclemass -
%   * mass of muscle, mass = PCSA*rho*lMopt, PCSA = Fmax/sigma
%       Note: rho is density (1058.7 kg/m³) and sigma is tension
%
%   - pctst -
%   * percentage of slow twitch fibers (0-1)
%
%   - vcemax -
%   * maximal muscle fiber velocities (default is 10*lMopt)
%
%   - Fiso -
%   * normalized muscle forces from active f-l relationship (FMltilde)
%
%   - Fmax -
%   * maximal isometric muscle forces
%
%   - modelmass -
%   * mass of the musculoskeletal model
%
%   - b -
%   * parameter determining transition smoothness for tanh approximations
%
%   - strength -
%   * muscle strength scaling factor
%
% OUTPUT:
%   - energy_total -
%   * total metabolic energy rate
%       energy_total: energy rate from totalHeatRate and Wdot
%       totalHeatRate: energy rate from Adot, Mdot and Sdot
%       Note: energy_total might be different than the sum of the different
%             components if the total heat rate was clamped to one.
%   
%   - Adot -
%   * energy rate from activation
%   
%   - Mdot -
%   * energy rate from maintenance
%   
%   - Sdot -
%   * energy rate from shortening and lengthening
%   
%   - Wdot -
%   * energy rate from positive mechanical work 
%   
%   - energy_model -
%   * energy rate from energy_total including basal rate
%
% Original author: Antoine Falisse
% Original date: 01/22/2019
%
% --------------------------------------------------------------------------
% This file is part of PredSim.
% 
% PredSim: A Framework for Rapid Predictive Simulations of Locomotion
% Copyright (c) 2026 KU Leuven
% 
% PredSim is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Affero General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% 
% PredSim is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public 
% License for more details.
% 
% You should have received a copy of the GNU Affero General Public License 
% along with PredSim. If not, see <https://www.gnu.org/licenses/>.
% --------------------------------------------------------------------------

%% scale the muscle mass to strength
musclemass=musclemass.*strength;

%% Parameters
% Ratio of slow twitch (st) and fast twitch (ft) fibers
pctft = 1-pctst;      

%% Twitch excitations
st_e = pctst.*sin(pi/2*exc);
ft_e = pctft.*(1-cos(pi/2*exc));

%% Activation heat rate
decay_function_value = 1;
activation_constant_st = 40; % default (Bhargava et al. 2004)
activation_constant_ft = 133; % default (Bhargava et al. 2004)
Adot = musclemass.*decay_function_value.*((activation_constant_st*st_e)+...
    (activation_constant_ft*ft_e));

%% Maintenance heat rate
fiber_length_dep = lMtilde;
maintenance_constant_st = 74; % default (Bhargava et al. 2004)
maintenance_constant_ft = 111; % default (Bhargava et al. 2004)
Mdot = musclemass.*fiber_length_dep.*((maintenance_constant_st*st_e)+...
    (maintenance_constant_ft*ft_e));

%% Shortening heat rate
fiber_force_total = Fce + Fpass;

alpha = 0.25*fiber_force_total;

% vM is positive (muscle lengthening)
vM_pos = 0.5 + 0.5*tanh(b*(vM));

% vM is negative (muscle shortening)

alpha = alpha + (-alpha + 0.157*fiber_force_total).*vM_pos;

Sdot = -alpha.*vM;

%% Mechanical work rate
% No negative mechanical work allowed
Wdot = - Fce.*vM;

%% Simple sum
totalHeatRate = Adot + Mdot + Sdot;

%% Total metabolic energy rate
energy_total = totalHeatRate + Wdot;

%% Energy model
basal_coef = 1.2; % default in OpenSim
basal_exp = 1; % default in OpenSim
energy_model = basal_coef*modelmass^basal_exp + sum(energy_total); 

end