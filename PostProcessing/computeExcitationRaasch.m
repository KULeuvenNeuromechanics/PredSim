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
% Last edit by: 
% Last edit date: 
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
