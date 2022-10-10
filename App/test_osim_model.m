
clear
clc

import org.opensim.modeling.*;

[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
osimpath = fullfile(pathRepo,'Subjects','Vitruvian_Man','Vitruvian_Man_wMarkers.osim');

model = Model(osimpath);
s = model.initSystem;


pos_head = model.getMarkerSet().get('top_of_head').getLocationInGround(s).getAsMat;
pos_shoulder = model.getBodySet().get('humerus_r').findBaseFrame().getPositionInGround(s).getAsMat;
pos_elbow = model.getBodySet().get('ulna_r').findBaseFrame().getPositionInGround(s).getAsMat;
pos_fingertip = model.getMarkerSet().get('fingertip_r').getLocationInGround(s).getAsMat;
pos_hip = model.getBodySet().get('femur_r').findBaseFrame().getPositionInGround(s).getAsMat;
pos_knee = model.getBodySet().get('tibia_r').findBaseFrame().getPositionInGround(s).getAsMat;
pos_sole = model.getMarkerSet().get('sole_r').getLocationInGround(s).getAsMat;
pos_toe = model.getMarkerSet().get('toe_r').getLocationInGround(s).getAsMat;

pos_knee(2) = pos_knee(2)-0.03;

%%

dim_height = pos_head(2)-pos_sole(2);
dim_fingertip_elbow = norm(pos_elbow-pos_fingertip);
dim_elbow_shoulder = norm(pos_elbow-pos_shoulder);
dim_shoulder_width = pos_shoulder(3)*2;
dim_hip_knee = pos_hip(2)-pos_knee(2);
dim_knee_ground = pos_knee(2)-pos_sole(2);
dim_foot_length = pos_toe(1)-pos_sole(1);


%%
dim_ratio_fingertip_elbow = dim_fingertip_elbow/dim_height;
dim_ratio_elbow_shoulder = dim_elbow_shoulder/dim_height;
dim_ratio_shoulder_width = dim_shoulder_width/dim_height;
dim_ratio_hip_knee = dim_hip_knee/dim_height;
dim_ratio_knee_ground = dim_knee_ground/dim_height;
dim_ratio_foot_length = dim_foot_length/dim_height;

(dim_ratio_fingertip_elbow + dim_ratio_elbow_shoulder)*2 + dim_ratio_shoulder_width
dim_ratio_knee_ground + dim_ratio_hip_knee

%%

sf_fingertip_elbow = (1/4*1.8)/dim_fingertip_elbow;
sf_elbow_shoulder = (1/8*1.8)/dim_elbow_shoulder;
sf_shoulder_width = (1/4*1.8)/dim_shoulder_width;
sf_hip_knee = (1/4*1.8)/dim_hip_knee;
sf_knee_ground = (1/4*1.8)/dim_knee_ground;
sf_foot_length = (1/7*1.8)/dim_foot_length;



