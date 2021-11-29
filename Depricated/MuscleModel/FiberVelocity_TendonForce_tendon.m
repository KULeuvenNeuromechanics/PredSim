% This function computes muscle fiber velocities from muscle-tendon forces.
% More details in De Groote et al. (2016): DOI: 10.1007/s10439-016-1591-9
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function [vM,vMtilde] = FiberVelocity_TendonForce_tendon(FTtile,...
    dFTtilde,params,lMT,vMT,Atendon,shift)

lMo = ones(size(FTtile,1),1)*params(2,:);
lTs = ones(size(FTtile,1),1)*params(3,:);
alphao = ones(size(FTtile,1),1)*params(4,:);
vMmax = ones(size(FTtile,1),1)*params(5,:);

% Inverse tendon force-length characteristic
lTtilde = log(5*(FTtile + 0.25 - shift))./Atendon + 0.995;

% Hill-type muscle model: geometric relationships
lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilde).^2);

vT = lTs.*dFTtilde./(0.2*Atendon*exp(Atendon*(lTtilde-0.995)));
cos_alpha = (lMT-lTs.*lTtilde)./lM;
vM = (vMT-vT).*cos_alpha;
vMtilde = vM./vMmax;

end
