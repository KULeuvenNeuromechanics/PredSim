% This function computes muscle fiber lengths from muscle-tendon forces.
% More details in De Groote et al. (2016): DOI: 10.1007/s10439-016-1591-9
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function [lM,lMtilde ] = ...
    FiberLength_TendonForce_tendon(FTtilde,params,lMT,aTendon,shift)

lMo = ones(size(FTtilde,1),1)*params(2,:);
lTs = ones(size(FTtilde,1),1)*params(3,:);
alphao = ones(size(FTtilde,1),1)*params(4,:);

% Tendon force-length characteristic
lTtilde = (log(5*(FTtilde + 0.25 - shift))/aTendon + 0.995);

% Hill-type muscle model: geometric relationships
lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilde).^2);
lMtilde = lM./lMo;

end
