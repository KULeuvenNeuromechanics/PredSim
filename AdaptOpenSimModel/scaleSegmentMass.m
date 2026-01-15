% --------------------------------------------------------------------------
% scaleSegmentMass
%   Scale the segment masses of an OpenSim model to get the desired mass 
%   distribution. 
% 
%
% Original author: Lars D'Hondt
% Original date: 15 May 2024
% --------------------------------------------------------------------------

clear
close all
clc


% original model
modelPath = "C:\GBW_MyPrograms\PredSim-dev\Subjects\DHondt_et_al_2024_3seg\DHondt_et_al_2024_3seg.osim";

% new model, with scaled mass
modelNewPath = "C:\GBW_MyPrograms\PredSim-dev\Subjects\DHondt_et_al_2025\DHondt_et_al_2025_v2.2.osim";


% For each segment, set the names of the bodies that are included and the
% desired percentage of body mass.

% HAT
segments(1).bodyNames = {'torso','pelvis'};
segments(1).massPct = 50.49;
% upper arm
segments(end+1).bodyNames = {'humerus_r','humerus_l'};
segments(end).massPct = 5.3167;
% lower arm
segments(end+1).bodyNames = {'radius_r','radius_l','ulna_r','ulna_l'};
segments(end).massPct = 3.1782;
% hand
segments(end+1).bodyNames = {'hand_r','hand_l'};
segments(end).massPct = 1.1967;
% upper leg
segments(end+1).bodyNames = {'femur_r','femur_l'};
segments(end).massPct = 24.3308;
% lower leg
segments(end+1).bodyNames = {'tibia_r','tibia_l'};
segments(end).massPct = 9.6956;
% foot
segments(end+1).bodyNames = {'talus_r','talus_l','calcn_r','calcn_l',...
    'midfoot_r','midfoot_l','forefoot_r','forefoot_l','toes_r','toes_l'};
segments(end).massPct = 5.7919;


%%
[mtot] = getModelMass(modelPath);

import org.opensim.modeling.*

model = Model(modelPath);

massVec = zeros(2,length(segments));

for j=1:length(segments)
    
    segMass = 0;
    for i=1:length(segments(j).bodyNames)
        body_i = model.getBodySet().get(segments(j).bodyNames{i});
        segMass = segMass + body_i.getMass();
    end

    sf_j = (mtot*segments(j).massPct/100) / segMass *62/mtot;
    scaleVec_j = Vec3(sf_j);
    
    for i=1:length(segments(j).bodyNames)
        body_i = model.getBodySet().get(segments(j).bodyNames{i});
        massVec(1,j) = massVec(1,j) + body_i.getMass();

        body_i.scaleMass(sf_j);
        massVec(2,j) = massVec(2,j) + body_i.getMass();

    end
end

massTotal = sum(massVec,2);
massPct = massVec./massTotal*100;


model.initSystem();
model.print(modelNewPath);



