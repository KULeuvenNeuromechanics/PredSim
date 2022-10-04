function fixContactSpherePositionAfterScaling(genericModelPath,scaledModelPath)
% OpenSim does not scale the position of the contact spheres during the
% scaling step. This code uses the change in distance between toes and
% calcaneus origin to scale the location of the contact spheres. Run this
% code after the scaling step.
% 
% INPUT:
%   - genericModelPath: Path to the generic model you scaled
%
%   - scaledModelPath: Path to the scaled model
% 
% OUTPUT:
%   - No output, the scaled model is overwritten with the corrected
%   locations of the contact spheres
% 
% Original author: Dhruv Gupta
% Original date: 04/October/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

%% Step 1: Get the position of toes origin with respect to calcaneus
import org.opensim.modeling.*;
model_1 = Model(genericModelPath);
toesInCalcn_r_1 = model_1.getJointSet.get('mtp_r').get_frames(0).get_translation.getAsMat';
toesInCalcn_l_1 = model_1.getJointSet.get('mtp_l').get_frames(0).get_translation.getAsMat';

model_2 = Model(scaledModelPath);
toesInCalcn_r_2 = model_2.getJointSet.get('mtp_r').get_frames(0).get_translation.getAsMat';
toesInCalcn_l_2 = model_2.getJointSet.get('mtp_l').get_frames(0).get_translation.getAsMat';

factorR = sqrt(sum(toesInCalcn_r_2.^2))/sqrt(sum(toesInCalcn_r_1.^2));
factorL = sqrt(sum(toesInCalcn_l_2.^2))/sqrt(sum(toesInCalcn_l_1.^2));
factorRX = toesInCalcn_r_2(1)/toesInCalcn_r_1(1);
factorLX = toesInCalcn_l_2(1)/toesInCalcn_l_1(1);

%% Step 2: Get the location of spheres in the generic model
clear cGS
cGS = model_1.get_ContactGeometrySet;
for i=1:cGS.getSize-1
    clear sphere
    sphere = cGS.get(i);
    sphereName{i} = char(sphere.getName);
    clear location
    location = sphere.get_location;
    sphereLocation(i,:) = [location.get(0) location.get(1) location.get(2)];
    if strfind(sphereName{i},'_r')
        scaledLocation(i,:) = factorR*sphereLocation(i,:);
    else
        scaledLocation(i,:) = factorL*sphereLocation(i,:);
    end        
end

%% Step 3: Put the scaled sphere locations in the scaled model and overwrite the model
clear cGS
cGS = model_2.get_ContactGeometrySet;
for i=1:cGS.getSize-1
    clear sphere
    sphere = cGS.get(i);
    idx = find(ismember(sphereName,char(sphere.getName)));
    clear location
    sphere.set_location(Vec3(scaledLocation(idx,1),scaledLocation(idx,2),scaledLocation(idx,3)));
end

model_2.initSystem();
model_2.print(scaledModelPath);
