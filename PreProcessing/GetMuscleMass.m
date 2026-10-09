function [massM,specific_tension] = GetMuscleMass(FMo,lMo,specific_tension)
% --------------------------------------------------------------------------
% GetMuscleMass
%   This function computes the mass of a muscle.
%   
% INPUT:
%   - FMo -
%   * Optimal active muscle force
%
%   - lMo -
%   * Optimal muscle fiber length
% 
%   - specific_tension -
%   * Specific tension of muscle fibers
%
% OUTPUT:
%   - massM -
%   * mass of a muscle
%
%   - specific_tension -
%   * Specific tension of muscle fibers
% 
% Original author: Dhruv Gupta
% Original date: 
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

volM = FMo.*lMo;
massM = volM.*(1059.7)./(specific_tension*1e6);

end

