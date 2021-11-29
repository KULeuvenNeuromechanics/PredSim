function [Names,params,lOpt,L_TendonSlack,Fiso,PennationAngle]=DispMusclesOsimModel(ModelPath)
% input= path to model and cell array with muscle names
% output= params (5xNMuscles) with  row: (1)  IsomForce (2)OptFiberLength
% 			(3) TendonSlackLength (4) PennationAngle (5) MaxFiberVelocity

% read the model
import org.opensim.modeling.*;
model = Model(ModelPath);

muscles = model.getMuscles();

for i = 1:muscles.getSize
   muscle = muscles.get(i-1);
   params(3,i) = muscle.getTendonSlackLength();		
   params(2,i) = muscle.getOptimalFiberLength(); 	
   params(1,i) = muscle.getMaxIsometricForce();  	
   params(4,i) = muscle.getPennationAngleAtOptimalFiberLength(); 
   params(5,i) = muscle.getMaxContractionVelocity()*params(2,i);
   Names{i} = char(muscle.getName());
end


% create additional variables with the same information
Fiso=params(1,:);
lOpt=params(2,:);
L_TendonSlack=params(3,:);
PennationAngle=params(4,:);

% print params
for i = 1:muscles.getSize
    disp([num2str(i) ' ' Names{i} ':     Fiso ' num2str(Fiso(i)) '      lOpt ' num2str(lOpt(i)) '     lTs ' num2str(L_TendonSlack(i)) ]);
end

end