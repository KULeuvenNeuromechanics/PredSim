function [energy_total,Adot,Mdot,Sdot,Wdot,energy_model] = ...
    getMetabolicEnergySmooth2004all_inclStrength(exc,act,lMtilde,vM,Fce,Fpass,...
        musclemass,pctst,Fiso,Fmax,modelmass,b,strength)
% --------------------------------------------------------------------------
% getMetabolicEnergySmooth2004all
%    This function computes the muscle energy expenditure based on the model
%    of Bhargava et al. (2004). This function provides a smooth approximation of
%    the original functions using tanh approximations for use in optimization
%    problems.
%
%    Note: we did not include the length dependence of the maintenance heat rate.
%      In OpenSim, they use a Piecewise Linear Function
%      https://github.com/opensim-org/opensim-core/blob/master/OpenSim/Simulation/Model/Bhargava2004MuscleMetabolicsProbe.cpp#L103
%    Note: we consider some lengthening heat rate, i.e. 
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
% Last edit by: 
% Last edit date: 
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
% F_iso that would be developed at the current activation and fiber length
% under isometric conditions (different than models of Umberger and Uchida).
% Fiso is getActiveForceLengthMultiplier in OpenSim. To minimize the
% difference between the models, we keep the same input, i.e., Fiso, that is
% used in the models of Umberger and Uchida.
F_iso = act.*Fiso.*Fmax.*strength; 
fiber_force_total = Fce + Fpass;
alpha = (0.16 * F_iso) + (0.18 * fiber_force_total);
% vM is positive (muscle lengthening)
vM_pos = 0.5 + 0.5*tanh(b*(vM));
% vM is negative (muscle shortening)
vM_neg = 1-vM_pos;
alpha = alpha + (-alpha + 0.157*fiber_force_total).*vM_pos;
% If get_use_force_dependent_shortening_prop_constant is false then
% alpha = 0.25*fiber_force_total;
% alpha(vM>0) = 0;
Sdot = -alpha.*vM;

%% Mechanical work rate
% No negative mechanical work allowed
Wdot = - Fce.*vM.*vM_neg;

%% Total power is non-negative
% If necessary, increase the shortening heat rate
% https://github.com/opensim-org/opensim-core/blob/master/OpenSim/Simulation/Model/Bhargava2004MuscleMetabolicsProbe.cpp#L393
Edot_W_beforeClamp = Adot + Mdot + Sdot + Wdot;
Edot_Wkg_beforeClamp_neg = 0.5+(0.5*tanh(b*(-Edot_W_beforeClamp)));
Sdot = Sdot - Edot_W_beforeClamp.*Edot_Wkg_beforeClamp_neg;

%% The total heat rate is not allowed to drop below 1(W/kg)
% https://github.com/opensim-org/opensim-core/blob/master/OpenSim/Simulation/Model/Bhargava2004MuscleMetabolicsProbe.cpp#L400
totalHeatRate = Adot + Mdot + Sdot;
% We first express in W/kg 
totalHeatRate = totalHeatRate./musclemass;
totalHeatRate = totalHeatRate + (-totalHeatRate + ...
    ones(size(totalHeatRate,1),1)).*(0.5 + 0.5*tanh(b*(1-totalHeatRate)));
% We then express back in W
totalHeatRate = totalHeatRate.*musclemass;

%% Total metabolic energy rate
energy_total = totalHeatRate + Wdot;

%% Energy model
basal_coef = 1.2; % default in OpenSim
basal_exp = 1; % default in OpenSim
energy_model = basal_coef*modelmass^basal_exp + sum(energy_total); 

end