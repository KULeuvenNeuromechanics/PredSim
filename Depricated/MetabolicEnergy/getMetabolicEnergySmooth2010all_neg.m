% This function computes the muscle energy expenditure based on the model
% of Umberger (2010). This function provides a smooth approximation of
% the original functions using tanh approximations for use in optimization
% problems. 
% This version treats negative mechanical work as Umberger et al. (2003).
%
% Author: Antoine Falisse
% Date: 01/22/2019
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
%   vcemax: maximal muscle fiber velocities (default is 10*lMopt)
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

function [energy_total,energy_am,energy_sl,energy_mech,energy_model] = ...
    getMetabolicEnergySmooth2010all_neg(exc,act,lMtilde,vMtilde,vM,Fce, ...
        musclemass,pctst,vcemax,Fiso,modelmass,b)

%% Parameters
% Percentage of slow twitch (st) and fast twitch (ft) fibers
pctft = 1-pctst;
% Scaling factor: aerobic (1.5) or anaerobic (1) activities
s = 1.5;            

%% Excitations and activations
a = exc + (-exc+act)/2.*(0.5+0.5*tanh(b*(act-exc)));

%% Heat activation and maintenance
% lMtilde is on the descending leg of the active force-length relationship
% i.e. lMtilde > 1 || lM > lMopt
lMtilde_des = 0.5 + 0.5*tanh(b*(lMtilde-1*ones(size(lMtilde,1),1)));        
hdotam = 128*pctft+25 + (-0.6*(128*pctft+25) + 0.6*(128*pctft+25) ... 
         .*Fiso).*lMtilde_des;       

aam         = a.^(0.6);
energy_am   = hdotam.*aam*s;

%% Heat shortening and lengthening
% vM is positive (muscle lengthening)
vMtilde_pos = 0.5 + 0.5*tanh(b*(vMtilde));
% vM is negative (muscle shortening)
vMtilde_neg = 1-vMtilde_pos;

coef_hs_ft  = 1*153.*ones(size(vcemax,1),1)./vcemax;
coef_hs_st  = 4*25.*ones(size(vcemax,1),1)./(vcemax/2.5); 
coef_hl     = 0.3*coef_hs_st; % different as compared to Umberger et al. (2003)

hdotsl = coef_hl.*vMtilde + (-coef_hl.*vMtilde ...
    - coef_hs_st.*vMtilde.*pctst ...
    - coef_hs_ft.*vMtilde.*pctft).*vMtilde_neg;

energy_sl = (hdotsl.*a.*s + ...
    (-hdotsl.*a.*s + hdotsl.*(a.^2)*s).*vMtilde_neg).*(1-lMtilde_des) + ...
    (hdotsl.*a.*s.*Fiso + (-hdotsl.*a.*s.*Fiso + ...
    hdotsl.*(a.^2)*s.*Fiso).*vMtilde_neg).*(lMtilde_des);

%% Mechanical work: same as Umberger et al. (2003)
% Include negative mechanical work
energy_mech = -1*Fce.*vM./musclemass;

%% The total rate of energy liberation is not allowed to drop below 1(W/kg)
% https://github.com/opensim-org/opensim-core/blob/master/OpenSim/Simulation/Model/Umberger2010MuscleMetabolicsProbe.cpp#L432
totalHeatRate = energy_am + energy_sl;
totalHeatRate = totalHeatRate + (-totalHeatRate + ...
    ones(size(totalHeatRate,1),1)).*(0.5 + 0.5*tanh(b*(1-totalHeatRate)));

%% Account for muscle mass
energy_am = energy_am.*musclemass;
energy_sl = energy_sl.*musclemass;
energy_mech = energy_mech.*musclemass;
energy_total = totalHeatRate.*musclemass+energy_mech;

%% Energy model
basal_coef = 1.2; % default in OpenSim
basal_exp = 1; % default in OpenSim
energy_model = basal_coef*modelmass^basal_exp + sum(energy_total);

end
