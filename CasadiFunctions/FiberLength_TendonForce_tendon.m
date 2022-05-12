function [lM,lMtilde,varargout] = ...
    FiberLength_TendonForce_tendon(FTtilde,lMo_in,lTs_in,alphao_in,lMT,aTendon,shift,MuscMoAsmp)
% --------------------------------------------------------------------------
% FiberLength_TendonForce_tendon
%    This function computes muscle fiber lengths from muscle-tendon forces.
%    More details in De Groote et al. (2016): DOI: 10.1007/s10439-016-1591-9
%   
% INPUT:
%
% OUTPUT:
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

lMo = ones(size(FTtilde,1),1)*lMo_in;
lTs = ones(size(FTtilde,1),1)*lTs_in;
alphao = ones(size(FTtilde,1),1)*alphao_in;

% Tendon force-length characteristic
lTtilde = (log(5*(FTtilde + 0.25 - shift))/aTendon + 0.995);

% Hill-type muscle model: geometric relationships
if(MuscMoAsmp == 0) % b = cst
    lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilde).^2);
else    % alpha = cst = alphao
    lM = (lMT-lTs.*lTtilde)./cos(alphao);
end
lMtilde = lM./lMo;
if nargout == 3
    varargout{1} = lTs.*lTtilde;
end
end