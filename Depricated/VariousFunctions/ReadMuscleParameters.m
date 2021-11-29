function [params,lOpt,L_TendonSlack,Fiso,PennationAngle]=ReadMuscleParameters(ModelPath,names)
% input= path to model and cell array with muscle names
% output= params (5xNMuscles) with  row: (1)  IsomForce (2)OptFiberLength
% 			(3) TendonSlackLength (4) PennationAngle (5) MaxFiberVelocity

% read the model
import org.opensim.modeling.*;
model = Model(ModelPath);

% read the muscle properties
nNames = length(names);
params = zeros(5, nNames);
muscles = model.getMuscles();

for i = 1:nNames
   muscle = muscles.get(names{i});
   params(3,i) = muscle.getTendonSlackLength();		
   params(2,i) = muscle.getOptimalFiberLength(); 	
   params(1,i) = muscle.getMaxIsometricForce();  	
   params(4,i) = muscle.getPennationAngleAtOptimalFiberLength(); 
   params(5,i) = muscle.getMaxContractionVelocity()*params(2,i);
end


% create additional variables with the same information
Fiso=params(1,:);
lOpt=params(2,:);
L_TendonSlack=params(3,:);
PennationAngle=params(4,:);
end