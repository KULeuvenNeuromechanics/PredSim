function [vM,vMtilde,varargout] = FiberVelocity_TendonForce_tendon(FTtile,...
    dFTtilde,lMo_in,lTs_in,alphao_in,vMmax_in,lMT,vMT,Atendon,shift,MuscMoAsmp)
% --------------------------------------------------------------------------
% FiberVelocity_TendonForce_tendon
%    This function computes muscle fiber velocities from muscle-tendon forces.
%    More details in De Groote et al. (2016): DOI: 10.1007/s10439-016-1591-9
%   
% INPUT:
%
% OUTPUT:
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
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

lMo = ones(size(FTtile,1),1)*lMo_in;
lTs = ones(size(FTtile,1),1)*lTs_in;
alphao = ones(size(FTtile,1),1)*alphao_in;
vMmax = ones(size(FTtile,1),1)*vMmax_in;

% Inverse tendon force-length characteristic
lTtilde = log(5*(FTtile + 0.25 - shift))./Atendon + 0.995;

% Hill-type muscle model: geometric relationships
vT = lTs.*dFTtilde./(0.2*Atendon*exp(Atendon*(lTtilde-0.995)));
if(MuscMoAsmp == 0) % b = cst
    lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilde).^2);
    cos_alpha = (lMT-lTs.*lTtilde)./lM;
else    % alpha = cst = alphao
    cos_alpha = cos(alphao);
end
vM = (vMT-vT).*cos_alpha;
vMtilde = vM./vMmax;

if nargout == 3
    varargout{1} = vT;
end

end