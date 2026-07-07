function [mtot] = getModelMass(ModelFile)
% --------------------------------------------------------------------------
% getModelMass
%   Computes to sum of all segment masses of an opensim model. 
%   Copied from https://github.com/KULeuvenNeuromechanics/NeuromechanicsToolkit/blob/75af1dead6960936d323b862b463ff789fb08692/OpenSimAPI/getModelMass.m
% 
% INPUT:
%   - Modelfile -
%   * path to the opensim model
%
% OUTPUT:
%   - mtot -
%   * total mass of the opensim model
% 
% Original author: Maarten Afschrift
% Original date: 6/Dec/2021
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

import org.opensim.modeling.*
m = Model(ModelFile);
mtot = 0;

for i=1:m.getBodySet.getSize()
    mtot = mtot + m.getBodySet.get(i-1).getMass();
end




end