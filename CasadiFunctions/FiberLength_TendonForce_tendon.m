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

lMo = ones(size(FTtilde,1),1)*lMo_in;
lTs = ones(size(FTtilde,1),1)*lTs_in;
alphao = ones(size(FTtilde,1),1)*alphao_in;

% Tendon force-length characteristic
lTtilde = (log(5*(FTtilde + 0.25 - shift))/aTendon + 0.995);

% Hill-type muscle model: geometric relationships
if(MuscMoAsmp == 0) % constantmuscle height/width
    lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilde).^2);
else    % constant pennation angle (= alphao)
    lM = (lMT-lTs.*lTtilde)./cos(alphao);
end
lMtilde = lM./lMo;

% Tendon length
if nargout == 3
    varargout{1} = lTs.*lTtilde;
end
end