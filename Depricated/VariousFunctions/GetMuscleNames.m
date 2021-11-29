function [Names]=GetMuscleNames(ModelPath)
% input= path to model and cell array with muscle names
% output= params (5xNMuscles) with  row: (1)  IsomForce (2)OptFiberLength
% 			(3) TendonSlackLength (4) PennationAngle (5) MaxFiberVelocity

% read the model
import org.opensim.modeling.*;
model = Model(ModelPath);

muscles = model.getMuscles();
Names = cell(1,muscles.getSize);
for i = 1:muscles.getSize
   muscle = muscles.get(i-1);
   Names{i} = char(muscle.getName());
end


end