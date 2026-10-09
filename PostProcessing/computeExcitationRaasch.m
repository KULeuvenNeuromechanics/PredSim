function e = computeExcitationRaasch(a, vA, tauDeact, tauAct)
% --------------------------------------------------------------------------
% computeExcitationRaasch
%   This function computes muscle excitations from time derivative of muscle 
%   activations using Raasch's model.
%   More details in De Groote et al. (2009): DOI: 10.1080/10255840902788587
% 
% INPUT:
%   - a -
%   * muscle activations
%
%   - vA -
%   * time derivatives of muscle activations
% 
%   - tauDeact -
%   * time constant of deactivation dynamics
%
%   - tauActct -
%   * time constant of activation dynamics
%
% OUTPUT:
%   - e -
%   * muscle excitations
% 
% Original author: Antoine Falisse
% Original date: 19/Dec/2018
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

td = (ones(size(a,1),1)*tauDeact);
ta = (ones(size(a,1),1)*tauAct);

e = zeros(size(a));
e(vA<=0) = td(vA<=0) .* vA(vA<=0) + a(vA<=0);

c1 = 1./ta - 1./td;
c2 = 1./td;
D = (c2 + c1 .* a).^2 + 4*c1.*vA;
e(vA>0) = (a(vA>0) .* c1(vA>0) - c2(vA>0) + sqrt(D(vA>0)))./(2*c1(vA>0));

end
