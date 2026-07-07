function [FMo, lMo, lTs, alphao, vMmax] = getMTparameters(modelPath,muscleNames)
% --------------------------------------------------------------------------
% getMTparameters
%   This function returns the muscle-tendon parameters of the muscles
%   specified in muscleNames from the model in modelPath.
%   
% INPUT:
%   - modelPath -
%   * path to the OpenSim model file (.osim)
% 
%   - muscleNames -
%   * muscles for which we want to read the muscle-tendon parameters
%
% OUTPUT:
%   - FMo -
%   * maximum isometric force (N)
%
%   - lMo -
%   * optimal fiber length (m)
%
%   - lTs -
%   * tendon slack length (m)
%
%   - alphao -
%   * pennation angle at optimal fiber length (rad)
%
%   - vMmax -
%   * maximum contraltion velocity (m/s)
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

import org.opensim.modeling.*

model = Model(modelPath);

NMuscles = length(muscleNames);
FMo = zeros(1,NMuscles);
lMo = zeros(1,NMuscles);
lTs = zeros(1,NMuscles);
alphao = zeros(1,NMuscles);
vMmax = zeros(1,NMuscles);

muscles = model.getMuscles();

for i = 1:NMuscles
   muscle = muscles.get(muscleNames{i});
   FMo(i) = muscle.getMaxIsometricForce();
   lMo(i) = muscle.getOptimalFiberLength();
   lTs(i) = muscle.getTendonSlackLength();
   alphao(i) = muscle.getPennationAngleAtOptimalFiberLength();
   vMmax(i) = muscle.getMaxContractionVelocity()*lMo(i);
end

end
