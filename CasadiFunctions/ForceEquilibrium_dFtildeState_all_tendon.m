function [dfse, dfM] = ...
    ForceEquilibrium_dFtildeState_all_tendon(a,lMtilde,vM,lT,FMo_in,lMo_in,...
    lTs_in,vMmax_in,Ftparam,Fvparam,Fpparam,Faparam,stiffness_shift,...
    stiffness_scale,strength,tendon_stiff)
% --------------------------------------------------------------------------
% ForceEquilibrium_dFtildeState_all_tendon
%    This function computes the partial derivative of the muscle and tendon
%    force curves, to compute the tendon and muscle stiffness.
%    Based on: De Groote et al. (2016): DOI: 10.1007/s10439-016-1591-9
%   
% INPUT:
%
% OUTPUT:
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
%   Adapted to compute derivatives by Menthy Denayer
% Last edit by: Menthy Denayer
% Last edit date: 03/June/2026
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


% define muscle properties
FMo = ones(size(a,1),1)*FMo_in;
lMo = ones(size(a,1),1)*lMo_in;
lTs = ones(size(a,1),1)*lTs_in;
vMmax = ones(size(a,1),1)*vMmax_in;

% Inverse tendon force-length characteristic
lTtilde = lT./lTs;

% define muscle lengths
vMtilde = vM./vMmax;

% tendon force-length characteristic
c1 = Ftparam(2);
c2 = Ftparam(3);

% tendon_stiff = kt*tendon_stiff;

dfse = c1 * tendon_stiff * exp((lTtilde - c2).*tendon_stiff);
dfse = dfse .* FMo ./ lTs;

% Active muscle force-length characteristic
b11 = Faparam(1);
b21 = Faparam(2);
b31 = Faparam(3);
b41 = Faparam(4);
b12 = Faparam(5);
b22 = Faparam(6);
b32 = Faparam(7);
b42 = Faparam(8);
b13 = 0.1;
b23 = 1;
b33 = 0.5*sqrt(0.5);
b43 = 0;

num3 = lMtilde-b23;
den3 = b33+b43*lMtilde;

ddlM3 = -0.5 * 2 * (den3 .* 1 - b43 * num3) ./ (den3.^2) .* (num3./den3);
dFMtilde3 = b13 * ddlM3 .* exp(-0.5*num3.^2./den3.^2);

num1 = lMtilde-b21;
den1 = b31+b41*lMtilde;

ddlM1 = -0.5 * 2 * (den1 .* 1 - b41 * num1) ./ (den1.^2) .* (num1./den1);
dFMtilde1 = b11 * ddlM1 .* exp(-0.5*num1.^2./den1.^2);

num2 = lMtilde-b22;
den2 = b32+b42*lMtilde;

ddlM2 = -0.5 * 2 * (den2 .* 1 - b42 * num2) ./ (den2.^2) .* (num2./den2);
dFMtilde2 = b12 * ddlM2 .* exp(-0.5*num2.^2./den2.^2);

dFMltilde = dFMtilde1 + dFMtilde2 + dFMtilde3;

% Passive muscle force-length characteristic
e0 = 0.6;
kpe = 4;

dt5 = kpe / (e0/stiffness_scale) * exp(kpe * (lMtilde - stiffness_shift) / (e0/stiffness_scale));

% Passive muscle force
dFpetilde = dt5 / Fpparam(2);

% Muscle force-velocity characteristic
d1 = Fvparam(1);
d2 = Fvparam(2);
d3 = Fvparam(3);
d4 = Fvparam(4);

FMvtilde = d1 * log((d2*vMtilde + d3) + sqrt((d2*vMtilde + d3).^2 + 1)) + d4;

% Active muscle force
dFcetilde = strength.*a.*dFMltilde.*FMvtilde;

% Total active & passive muscle force
dfM = dFcetilde + dFpetilde;
dfM = dfM .* FMo ./ lMo;

end
