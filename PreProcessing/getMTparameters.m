% This function returns the muscle-tendon parameters of the muscles
% specified in muscleNames from the model in modelPath.
%
% These parameters are:
%   - FMo       maximum isometric force (N)
%   - lMo       optimal fiber length (m)
%   - lTs       tendon slack length (m)
%   - alphao    pennation angle at optimal fiber length (rad)
%   - vMmax     maximum contraltion velocity (m/s)
%
% Author: Antoine Falisse
% Date: 12/19/2018
% Adapted by Lars D'Hondt on 17/March/2022 to output each parameter
% separately.
%
function [FMo, lMo, lTs, alphao, vMmax] = getMTparameters(modelPath,muscleNames)

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
