function fixContactSpherePositionAfterScaling(genericModelPath,scaledModelPath,varargin)
% This code has two applications:
% Application 1:
% OpenSim does not scale the position of the contact spheres during the
% scaling step. This code uses the change in distance between toes and
% calcaneus origin to scale the location of the contact spheres. Run this
% code after the scaling step.
% Application 2:
% This code can also be used to change the location of contact spheres of
% any model to represent scaled locations based off of the specified
% generic model iff the two models have the same number of spheres with the
% same names.
% 
% INPUT:
%   - genericModelPath: Path to the generic model.
%
%   - scaledModelPath: Path to the scaled model (appliaction 1) or the
%   model whose contact sphere locations you want to fix (application 2).
% 
%   - Optional last input: A 1 x 3 matrix with the scaling factors of
%   x-y-z locations of the spheres. If this input is provided, these inputs
%   are used to scale the location of the spheres. If one of the inputs is
%   nan, then the scaling factor based on the change in distance between
%   toes and calcaneus origin is used to scale the location of that
%   direction. For example, if the last input is [nan 1 nan], it will use 1
%   as the scaling factor of the contact sphere location in y direction and
%   the scaling factor based on the change in distance between toes and
%   calcaneus origin for the x and z direction.
% 
% OUTPUT:
%   - No output, the scaled model is overwritten with the corrected
%   locations of the contact spheres
% 
% Original author: Dhruv Gupta
% Original date: 04/October/2022
%
% Last edit by: Ellis Van Can   
% Last edit date: 10/February/2026
% --------------------------------------------------------------------------

if nargin>2
    factorR = varargin{1};
    factorL = varargin{2};
end

%% Step 1: Get the position of toes origin with respect to calcaneus
import org.opensim.modeling.*;
model_1 = Model(genericModelPath);
toesInCalcn_r_1 = model_1.getJointSet.get('mtp_r').get_frames(0).get_translation.getAsMat';
toesInCalcn_l_1 = model_1.getJointSet.get('mtp_l').get_frames(0).get_translation.getAsMat';

model_2 = Model(scaledModelPath);
toesInCalcn_r_2 = model_2.getJointSet.get('mtp_r').get_frames(0).get_translation.getAsMat';
toesInCalcn_l_2 = model_2.getJointSet.get('mtp_l').get_frames(0).get_translation.getAsMat';

if ~exist('factorR')
    factorR = ones(1,3)*sqrt(sum(toesInCalcn_r_2.^2))/sqrt(sum(toesInCalcn_r_1.^2));
    factorL = ones(1,3)*sqrt(sum(toesInCalcn_l_2.^2))/sqrt(sum(toesInCalcn_l_1.^2));
else
    factorR(isnan(factorR)) = sqrt(sum(toesInCalcn_r_2.^2))/sqrt(sum(toesInCalcn_r_1.^2));
    factorL(isnan(factorL)) = sqrt(sum(toesInCalcn_l_2.^2))/sqrt(sum(toesInCalcn_l_1.^2));
end

%% Step 2: Get the location of spheres in the generic model
clear cGS
cGS = model_1.get_ContactGeometrySet;
sphereName = {};  sphereLocation = [];  scaledLocation = [];
for i = 0:cGS.getSize-1
    contactGeom = cGS.get(i);
    
    % Check if it's a ContactSphere (not ContactHalfSpace for floor)
    if strcmp(char(contactGeom.getConcreteClassName()), 'ContactSphere')
        sphereName{end+1} = char(contactGeom.getName);
        sphereLocation(end+1,:) = contactGeom.get_location.getAsMat';
        
        if contains(char(contactGeom.getName), '_r')
            scaledLocation(end+1,:) = factorR.*sphereLocation(end,:);
        else
            scaledLocation(end+1,:) = factorL.*sphereLocation(end,:);
        end
    end
end

%% Step 3: Put the scaled sphere locations in the scaled model and overwrite the model
clear cGS
cGS = model_2.get_ContactGeometrySet;

for i = 0:cGS.getSize-1
    contactGeom = cGS.get(i);
    
    % Check if it's a ContactSphere
    if strcmp(char(contactGeom.getConcreteClassName()), 'ContactSphere')
        geomName = char(contactGeom.getName);
        idx = ismember(sphereName, geomName);        
        if any(idx)
            contactGeom.set_location(Vec3.createFromMat(scaledLocation(idx,:)));
        end
    end
end

model_2.initSystem();
model_2.print(scaledModelPath);

