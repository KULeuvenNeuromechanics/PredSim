% This function computes muscle excitations from time derivative of muscle 
% activations using Raasch's model.
% More details in De Groote et al. (2009): DOI: 10.1080/10255840902788587
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function e = computeExcitationRaasch(a, vA, tauDeact, tauAct)

td = (ones(size(a,1),1)*tauDeact);
ta = (ones(size(a,1),1)*tauAct);

e = zeros(size(a));
e(vA<=0) = td(vA<=0) .* vA(vA<=0) + a(vA<=0);

c1 = 1./ta - 1./td;
c2 = 1./td;
D = (c2 + c1 .* a).^2 + 4*c1.*vA;
e(vA>0) = (a(vA>0) .* c1(vA>0) - c2(vA>0) + sqrt(D(vA>0)))./(2*c1(vA>0));

end
