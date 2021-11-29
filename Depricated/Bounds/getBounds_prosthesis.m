% This script provides bounds and scaling factors for the design variables.
% The bounds on the joint variables are informed by experimental data.
% The bounds on the remaining variables are fixed.
% The bounds are scaled such that the upper/lower bounds cannot be
% larger/smaller than 1/-1.
%
% Author: Antoine Falisse
% Date: 12/19/2018
%
function [bounds,scaling] = getBounds_prosthesis(Qs,NMuscle,nq,jointi)

%% Spline approximation of Qs to get Qdots and Qdotdots
Qs_spline.data = zeros(size(Qs.allfilt));
Qs_spline.data(:,1) = Qs.allfilt(:,1);
Qdots_spline.data = zeros(size(Qs.allfilt));
Qdots_spline.data(:,1) = Qs.allfilt(:,1);
Qdotdots_spline.data = zeros(size(Qs.allfilt));
Qdotdots_spline.data(:,1) = Qs.allfilt(:,1);
for i = 2:size(Qs.allfilt,2)
    Qs.datafiltspline(i) = spline(Qs.allfilt(:,1),Qs.allfilt(:,i));
    [Qs_spline.data(:,i),Qdots_spline.data(:,i),...
        Qdotdots_spline.data(:,i)] = ...
        SplineEval_ppuval(Qs.datafiltspline(i),Qs.allfilt(:,1),1);
end

%% Qs
% The extreme values are selected as upper/lower bounds, which are then
% further extended.
% Pelvis tilt
bounds.Qs.upper(jointi.pelvis.tilt) = max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'))));
bounds.Qs.lower(jointi.pelvis.tilt) = min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'))));
% Pelvis list
bounds.Qs.upper(jointi.pelvis.list) = max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_list'))));
bounds.Qs.lower(jointi.pelvis.list) = min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_list'))));
% Pelvis tilt
bounds.Qs.upper(jointi.pelvis.rot) = max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rotation'))));
bounds.Qs.lower(jointi.pelvis.rot) = min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rotation'))));
% Hip flexion
bounds.Qs.upper(jointi.hip_flex.l) = max(max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r')))...
    ,max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'))));
bounds.Qs.lower(jointi.hip_flex.l) = min(min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'))));
bounds.Qs.upper(jointi.hip_flex.r) = bounds.Qs.upper(jointi.hip_flex.l);
bounds.Qs.lower(jointi.hip_flex.r) = bounds.Qs.lower(jointi.hip_flex.l);
% Hip adduction
bounds.Qs.upper(jointi.hip_add.l) = max(max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_l'))));
bounds.Qs.lower(jointi.hip_add.l) = min(min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_l'))));
bounds.Qs.upper(jointi.hip_add.r) = bounds.Qs.upper(jointi.hip_add.l);
bounds.Qs.lower(jointi.hip_add.r) = bounds.Qs.lower(jointi.hip_add.l);
% Hip rotation
bounds.Qs.upper(jointi.hip_rot.l) = max(max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_l'))));
bounds.Qs.lower(jointi.hip_rot.l) = min(min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_l'))));
bounds.Qs.upper(jointi.hip_rot.r) = bounds.Qs.upper(jointi.hip_rot.l);
bounds.Qs.lower(jointi.hip_rot.r) = bounds.Qs.lower(jointi.hip_rot.l);
% Knee
bounds.Qs.upper(jointi.knee.l) = max(max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'))));
bounds.Qs.lower(jointi.knee.l) = min(min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'))));
bounds.Qs.upper(jointi.knee.r) = bounds.Qs.upper(jointi.knee.l);
bounds.Qs.lower(jointi.knee.r) = bounds.Qs.lower(jointi.knee.l);
% Ankle
bounds.Qs.upper(jointi.ankle.l) = max(max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'))));
bounds.Qs.lower(jointi.ankle.l) = min(min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'))));
bounds.Qs.upper(jointi.ankle.r) = bounds.Qs.upper(jointi.ankle.l);
bounds.Qs.lower(jointi.ankle.r) = bounds.Qs.lower(jointi.ankle.l);
% Subtalar
bounds.Qs.upper(jointi.subt.l) = max(max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_l'))));
bounds.Qs.lower(jointi.subt.l) = min(min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_l'))));
bounds.Qs.upper(jointi.subt.r) = bounds.Qs.upper(jointi.subt.l);
bounds.Qs.lower(jointi.subt.r) = bounds.Qs.lower(jointi.subt.l);
% Trunk extension
bounds.Qs.upper(jointi.trunk.ext) = max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'))));
bounds.Qs.lower(jointi.trunk.ext) = min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'))));
% Trunk bending
bounds.Qs.upper(jointi.trunk.ben) = max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_bending'))));
bounds.Qs.lower(jointi.trunk.ben) = min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_bending'))));
% Trunk rotation
bounds.Qs.upper(jointi.trunk.rot) = max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_rotation'))));
bounds.Qs.lower(jointi.trunk.rot) = min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_rotation'))));
% Shoulder flexion
bounds.Qs.upper(jointi.sh_flex.l) = max(max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_r')))...
    ,max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_l'))));
bounds.Qs.lower(jointi.sh_flex.l) = min(min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_l'))));
bounds.Qs.upper(jointi.sh_flex.r) = bounds.Qs.upper(jointi.sh_flex.l);
bounds.Qs.lower(jointi.sh_flex.r) = bounds.Qs.lower(jointi.sh_flex.l);
% Shoulder adduction
bounds.Qs.upper(jointi.sh_add.l) = max(max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_l'))));
bounds.Qs.lower(jointi.sh_add.l) = min(min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_l'))));
bounds.Qs.upper(jointi.sh_add.r) = bounds.Qs.upper(jointi.sh_add.l);
bounds.Qs.lower(jointi.sh_add.r) = bounds.Qs.lower(jointi.sh_add.l);
rec_uw_sh_add = bounds.Qs.upper(jointi.sh_add.l);
% Shoulder rotation
bounds.Qs.upper(jointi.sh_rot.l) = max(max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'))));
bounds.Qs.lower(jointi.sh_rot.l) = min(min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'))));
bounds.Qs.upper(jointi.sh_rot.r) = bounds.Qs.upper(jointi.sh_rot.l);
bounds.Qs.lower(jointi.sh_rot.r) = bounds.Qs.lower(jointi.sh_rot.l);
rec_uw_sh_rot = bounds.Qs.upper(jointi.sh_rot.l);
% Elbow
bounds.Qs.upper(jointi.elb.l) = max(max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'))));
bounds.Qs.lower(jointi.elb.l) = min(min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'))));
bounds.Qs.upper(jointi.elb.r) = bounds.Qs.upper(jointi.elb.l);
bounds.Qs.lower(jointi.elb.r) = bounds.Qs.lower(jointi.elb.l);
% The bounds are extended by twice the absolute difference between upper
% and lower bounds.
Qs_range = abs(bounds.Qs.upper - bounds.Qs.lower);
bounds.Qs.lower = bounds.Qs.lower - 2*Qs_range;
bounds.Qs.upper = bounds.Qs.upper + 2*Qs_range;
% For several joints, we manually adjust the bounds
% Pelvis_tx
bounds.Qs.upper(jointi.pelvis.tx) = 2;  
bounds.Qs.lower(jointi.pelvis.tx) = 0;
% Pelvis_ty
bounds.Qs.upper(jointi.pelvis.ty) = 1.1;  
bounds.Qs.lower(jointi.pelvis.ty) = 0.75;
% Pelvis_tz
bounds.Qs.upper(jointi.pelvis.tz) = 0.1;
bounds.Qs.lower(jointi.pelvis.tz) = -0.1;
% Elbow
bounds.Qs.lower(jointi.elb.l) = 0;
bounds.Qs.lower(jointi.elb.r) = 0;
% Shoulder adduction
bounds.Qs.upper(jointi.sh_add.l) = rec_uw_sh_add;
bounds.Qs.upper(jointi.sh_add.r) = rec_uw_sh_add;
% Shoulder rotation
bounds.Qs.upper(jointi.sh_rot.l) = rec_uw_sh_rot;
bounds.Qs.upper(jointi.sh_rot.r) = rec_uw_sh_rot;

%% Qdots
% The extreme values are selected as upper/lower bounds, which are then
% further extended.
% Pelvis tilt
bounds.Qdots.upper(jointi.pelvis.tilt) = max((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'))));
bounds.Qdots.lower(jointi.pelvis.tilt) = min((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'))));
% Pelvis list
bounds.Qdots.upper(jointi.pelvis.list) = max((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_list'))));
bounds.Qdots.lower(jointi.pelvis.list) = min((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_list'))));
% Pelvis rotation
bounds.Qdots.upper(jointi.pelvis.rot) = max((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rotation'))));
bounds.Qdots.lower(jointi.pelvis.rot) = min((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rotation'))));
% Pelvis_tx
bounds.Qdots.upper(jointi.pelvis.tx) = max((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx')))); 
bounds.Qdots.lower(jointi.pelvis.tx) = min((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'))));
% Pelvis_ty
bounds.Qdots.upper(jointi.pelvis.ty) = max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'))); 
bounds.Qdots.lower(jointi.pelvis.ty) = min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'))); 
% Pelvis_tz
bounds.Qdots.upper(jointi.pelvis.tz) = max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz'))); 
bounds.Qdots.lower(jointi.pelvis.tz) = min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz')));
% Hip flexion
bounds.Qdots.upper(jointi.hip_flex.l) = max(max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'))));
bounds.Qdots.lower(jointi.hip_flex.l) = min(min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'))));
bounds.Qdots.upper(jointi.hip_flex.r) = bounds.Qdots.upper(jointi.hip_flex.l);
bounds.Qdots.lower(jointi.hip_flex.r) = bounds.Qdots.lower(jointi.hip_flex.l);
% Hip adduction
bounds.Qdots.upper(jointi.hip_add.l) = max(max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_r'))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_l'))));
bounds.Qdots.lower(jointi.hip_add.l) = min(min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_l'))));
bounds.Qdots.upper(jointi.hip_add.r) = bounds.Qdots.upper(jointi.hip_add.l);
bounds.Qdots.lower(jointi.hip_add.r) = bounds.Qdots.lower(jointi.hip_add.l);
% Hip rotation
bounds.Qdots.upper(jointi.hip_rot.l) = max(max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_r'))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_l'))));
bounds.Qdots.lower(jointi.hip_rot.l) = min(min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_l'))));
bounds.Qdots.upper(jointi.hip_rot.r) = bounds.Qdots.upper(jointi.hip_rot.l);
bounds.Qdots.lower(jointi.hip_rot.r) = bounds.Qdots.lower(jointi.hip_rot.l);
% Knee
bounds.Qdots.upper(jointi.knee.l) = max(max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'))));
bounds.Qdots.lower(jointi.knee.l) = min(min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'))));
bounds.Qdots.upper(jointi.knee.r) = bounds.Qdots.upper(jointi.knee.l);
bounds.Qdots.lower(jointi.knee.r) = bounds.Qdots.lower(jointi.knee.l);
% Ankle
bounds.Qdots.upper(jointi.ankle.l) = max(max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'))));
bounds.Qdots.lower(jointi.ankle.l) = min(min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'))));
bounds.Qdots.upper(jointi.ankle.r) = bounds.Qdots.upper(jointi.ankle.l);
bounds.Qdots.lower(jointi.ankle.r) = bounds.Qdots.lower(jointi.ankle.l);
% Subtalar
bounds.Qdots.upper(jointi.subt.l) = max(max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_r'))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_l'))));
bounds.Qdots.lower(jointi.subt.l) = min(min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_l'))));
bounds.Qdots.upper(jointi.subt.r) = bounds.Qdots.upper(jointi.subt.l);
bounds.Qdots.lower(jointi.subt.r) = bounds.Qdots.lower(jointi.subt.l);
% Trunk extension
bounds.Qdots.upper(jointi.trunk.ext) = max((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'))));
bounds.Qdots.lower(jointi.trunk.ext) = min((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'))));
% Trunk bending
bounds.Qdots.upper(jointi.trunk.ben) = max((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_bending'))));
bounds.Qdots.lower(jointi.trunk.ben) = min((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_bending'))));
% Trunk rotation
bounds.Qdots.upper(jointi.trunk.rot) = max((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_rotation'))));
bounds.Qdots.lower(jointi.trunk.rot) = min((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_rotation'))));
% Shoulder flexion
bounds.Qdots.upper(jointi.sh_flex.l) = max(max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_r')))...
    ,max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_l'))));
bounds.Qdots.lower(jointi.sh_flex.l) = min(min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_l'))));
bounds.Qdots.upper(jointi.sh_flex.r) = bounds.Qdots.upper(jointi.sh_flex.l);
bounds.Qdots.lower(jointi.sh_flex.r) = bounds.Qdots.lower(jointi.sh_flex.l);
% Shoulder adduction
bounds.Qdots.upper(jointi.sh_add.l) = max(max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_r'))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_l'))));
bounds.Qdots.lower(jointi.sh_add.l) = min(min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_l'))));
bounds.Qdots.upper(jointi.sh_add.r) = bounds.Qdots.upper(jointi.sh_add.l);
bounds.Qdots.lower(jointi.sh_add.r) = bounds.Qdots.lower(jointi.sh_add.l);
% Shoulder rotation
bounds.Qdots.upper(jointi.sh_rot.l) = max(max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'))));
bounds.Qdots.lower(jointi.sh_rot.l) = min(min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'))));
bounds.Qdots.upper(jointi.sh_rot.r) = bounds.Qdots.upper(jointi.sh_rot.l);
bounds.Qdots.lower(jointi.sh_rot.r) = bounds.Qdots.lower(jointi.sh_rot.l);
% Elbow
bounds.Qdots.upper(jointi.elb.l) = max(max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'))));
bounds.Qdots.lower(jointi.elb.l) = min(min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'))));
bounds.Qdots.upper(jointi.elb.r) = bounds.Qdots.upper(jointi.elb.l);
bounds.Qdots.lower(jointi.elb.r) = bounds.Qdots.lower(jointi.elb.l);
% The bounds are extended by 3 times the absolute difference between upper
% and lower bounds.
Qdots_range = abs(bounds.Qdots.upper - bounds.Qdots.lower);
bounds.Qdots.lower = bounds.Qdots.lower - 3*Qdots_range;
bounds.Qdots.upper = bounds.Qdots.upper + 3*Qdots_range;

%% Qdotdots
% The extreme values are selected as upper/lower bounds, which are then
% further extended.
% Pelvis tilt
bounds.Qdotdots.upper(jointi.pelvis.tilt) = max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt')));
bounds.Qdotdots.lower(jointi.pelvis.tilt) = min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt')));
% Pelvis list
bounds.Qdotdots.upper(jointi.pelvis.list) = max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_list')));
bounds.Qdotdots.lower(jointi.pelvis.list) = min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_list')));
% Pelvis rotation
bounds.Qdotdots.upper(jointi.pelvis.rot) = max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rotation')));
bounds.Qdotdots.lower(jointi.pelvis.rot) = min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rotation')));
% Pelvis_tx
bounds.Qdotdots.upper(jointi.pelvis.tx) = max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'))); 
bounds.Qdotdots.lower(jointi.pelvis.tx) = min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx')));
% Pelvis_ty
bounds.Qdotdots.upper(jointi.pelvis.ty) = max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty')));
bounds.Qdotdots.lower(jointi.pelvis.ty) = min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty')));
% Pelvis_tz
bounds.Qdotdots.upper(jointi.pelvis.tz) = max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz')));
bounds.Qdotdots.lower(jointi.pelvis.tz) = min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz')));
% Hip flexion
bounds.Qdotdots.upper(jointi.hip_flex.l) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'))));
bounds.Qdotdots.lower(jointi.hip_flex.l) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'))));
bounds.Qdotdots.upper(jointi.hip_flex.r) = bounds.Qdotdots.upper(jointi.hip_flex.l);
bounds.Qdotdots.lower(jointi.hip_flex.r) = bounds.Qdotdots.lower(jointi.hip_flex.l);
% Hip adduction
bounds.Qdotdots.upper(jointi.hip_add.l) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_r'))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_l'))));
bounds.Qdotdots.lower(jointi.hip_add.l) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_l'))));
bounds.Qdotdots.upper(jointi.hip_add.r) = bounds.Qdotdots.upper(jointi.hip_add.l);
bounds.Qdotdots.lower(jointi.hip_add.r) = bounds.Qdotdots.lower(jointi.hip_add.l);
% Hip rotation
bounds.Qdotdots.upper(jointi.hip_rot.l) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_r'))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_l'))));
bounds.Qdotdots.lower(jointi.hip_rot.l) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_l'))));
bounds.Qdotdots.upper(jointi.hip_rot.r) = bounds.Qdotdots.upper(jointi.hip_rot.l);
bounds.Qdotdots.lower(jointi.hip_rot.r) = bounds.Qdotdots.lower(jointi.hip_rot.l);
% Knee
bounds.Qdotdots.upper(jointi.knee.l) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'))));
bounds.Qdotdots.lower(jointi.knee.l) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'))));
bounds.Qdotdots.upper(jointi.knee.r) = bounds.Qdotdots.upper(jointi.knee.l);
bounds.Qdotdots.lower(jointi.knee.r) = bounds.Qdotdots.lower(jointi.knee.l);
% Ankle
bounds.Qdotdots.upper(jointi.ankle.l) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'))));
bounds.Qdotdots.lower(jointi.ankle.l) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'))));
bounds.Qdotdots.upper(jointi.ankle.r) = bounds.Qdotdots.upper(jointi.ankle.l);
bounds.Qdotdots.lower(jointi.ankle.r) = bounds.Qdotdots.lower(jointi.ankle.l);
% Subtalar
bounds.Qdotdots.upper(jointi.subt.l) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_r'))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_l'))));
bounds.Qdotdots.lower(jointi.subt.l) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_l'))));
bounds.Qdotdots.upper(jointi.subt.r) = bounds.Qdotdots.upper(jointi.subt.l);
bounds.Qdotdots.lower(jointi.subt.r) = bounds.Qdotdots.lower(jointi.subt.l);
% Trunk extension
bounds.Qdotdots.upper(jointi.trunk.ext) = max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension')));
bounds.Qdotdots.lower(jointi.trunk.ext) = min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension')));
% Trunk bending
bounds.Qdotdots.upper(jointi.trunk.ben) = max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_bending')));
bounds.Qdotdots.lower(jointi.trunk.ben) = min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_bending')));
% Trunk rotation
bounds.Qdotdots.upper(jointi.trunk.rot) = max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_rotation')));
bounds.Qdotdots.lower(jointi.trunk.rot) = min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_rotation')));
% Shoulder flexion
bounds.Qdotdots.upper(jointi.sh_flex.l) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_r')))...
    ,max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_l'))));
bounds.Qdotdots.lower(jointi.sh_flex.l) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_l'))));
bounds.Qdotdots.upper(jointi.sh_flex.r) = bounds.Qdotdots.upper(jointi.sh_flex.l);
bounds.Qdotdots.lower(jointi.sh_flex.r) = bounds.Qdotdots.lower(jointi.sh_flex.l);
% Shoulder adduction
bounds.Qdotdots.upper(jointi.sh_add.l) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_r'))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_l'))));
bounds.Qdotdots.lower(jointi.sh_add.l) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_l'))));
bounds.Qdotdots.upper(jointi.sh_add.r) = bounds.Qdotdots.upper(jointi.sh_add.l);
bounds.Qdotdots.lower(jointi.sh_add.r) = bounds.Qdotdots.lower(jointi.sh_add.l);
% Shoulder rotation
bounds.Qdotdots.upper(jointi.sh_rot.l) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'))));
bounds.Qdotdots.lower(jointi.sh_rot.l) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'))));
bounds.Qdotdots.upper(jointi.sh_rot.r) = bounds.Qdotdots.upper(jointi.sh_rot.l);
bounds.Qdotdots.lower(jointi.sh_rot.r) = bounds.Qdotdots.lower(jointi.sh_rot.l);
% Elbow angle
bounds.Qdotdots.upper(jointi.elb.l) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'))));
bounds.Qdotdots.lower(jointi.elb.l) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'))));
bounds.Qdotdots.upper(jointi.elb.r) = bounds.Qdotdots.upper(jointi.elb.l);
bounds.Qdotdots.lower(jointi.elb.r) = bounds.Qdotdots.lower(jointi.elb.l);
% The bounds are extended by 3 times the absolute difference between upper
% and lower bounds.
Qdotdots_range = abs(bounds.Qdotdots.upper - bounds.Qdotdots.lower);
bounds.Qdotdots.lower = bounds.Qdotdots.lower - 3*Qdotdots_range;
bounds.Qdotdots.upper = bounds.Qdotdots.upper + 3*Qdotdots_range;

%% Muscle activations
bounds.a.lower = 0.05*ones(1,NMuscle);
bounds.a.upper = ones(1,NMuscle);

%% Muscle-tendon forces
bounds.FTtilde.lower = zeros(1,NMuscle);
bounds.FTtilde.upper = 5*ones(1,NMuscle);

%% Time derivative of muscle activations
tact = 0.015;
tdeact = 0.06;
bounds.vA.lower = (-1/100*ones(1,NMuscle))./(ones(1,NMuscle)*tdeact);
bounds.vA.upper = (1/100*ones(1,NMuscle))./(ones(1,NMuscle)*tact);

%% Time derivative of muscle-tendon forces
bounds.dFTtilde.lower = -1*ones(1,NMuscle);
bounds.dFTtilde.upper = 1*ones(1,NMuscle);

%% Arm activations
bounds.a_a.lower = -ones(1,nq.arms);
bounds.a_a.upper = ones(1,nq.arms);

%% Arm excitations
bounds.e_a.lower = -ones(1,nq.arms);
bounds.e_a.upper = ones(1,nq.arms);

%% Final time
bounds.tf.lower = 0.6;
bounds.tf.upper = 2;

%% Scaling
% Qs
scaling.Qs      = max(abs(bounds.Qs.lower),abs(bounds.Qs.upper));
bounds.Qs.lower = (bounds.Qs.lower)./scaling.Qs;
bounds.Qs.upper = (bounds.Qs.upper)./scaling.Qs;
% Qdots
scaling.Qdots      = max(abs(bounds.Qdots.lower),abs(bounds.Qdots.upper));
bounds.Qdots.lower = (bounds.Qdots.lower)./scaling.Qdots;
bounds.Qdots.upper = (bounds.Qdots.upper)./scaling.Qdots;
% Qs and Qdots are intertwined
bounds.QsQdots.lower = zeros(1,2*nq.all);
bounds.QsQdots.upper = zeros(1,2*nq.all);
bounds.QsQdots.lower(1,1:2:end) = bounds.Qs.lower;
bounds.QsQdots.upper(1,1:2:end) = bounds.Qs.upper;
bounds.QsQdots.lower(1,2:2:end) = bounds.Qdots.lower;
bounds.QsQdots.upper(1,2:2:end) = bounds.Qdots.upper;
scaling.QsQdots                 = zeros(1,2*nq.all);
scaling.QsQdots(1,1:2:end)      = scaling.Qs ;
scaling.QsQdots(1,2:2:end)      = scaling.Qdots ;
% Qdotdots
scaling.Qdotdots = max(abs(bounds.Qdotdots.lower),...
    abs(bounds.Qdotdots.upper));
bounds.Qdotdots.lower = (bounds.Qdotdots.lower)./scaling.Qdotdots;
bounds.Qdotdots.upper = (bounds.Qdotdots.upper)./scaling.Qdotdots;
bounds.Qdotdots.lower(isnan(bounds.Qdotdots.lower)) = 0;
bounds.Qdotdots.upper(isnan(bounds.Qdotdots.upper)) = 0;
% Arm torque actuators
% Fixed scaling factor
scaling.ArmTau = 150;
% Time derivative of muscle activations
% Fixed scaling factor
scaling.vA = 100;
% Muscle activations
scaling.a = 1;
% Arm activations
scaling.a_a = 1;
% Arm excitations
scaling.e_a = 1;
% Time derivative of muscle-tendon forces
% Fixed scaling factor
scaling.dFTtilde = 100;
% Muscle-tendon forces
scaling.FTtilde         = max(...
    abs(bounds.FTtilde.lower),abs(bounds.FTtilde.upper)); 
bounds.FTtilde.lower    = (bounds.FTtilde.lower)./scaling.FTtilde;
bounds.FTtilde.upper    = (bounds.FTtilde.upper)./scaling.FTtilde;

%% Hard bounds
% We impose the initial position of pelvis_tx to be 0
bounds.QsQdots_0.lower = bounds.QsQdots.lower;
bounds.QsQdots_0.upper = bounds.QsQdots.upper;
bounds.QsQdots_0.lower(2*jointi.pelvis.tx-1) = 0;
bounds.QsQdots_0.upper(2*jointi.pelvis.tx-1) = 0;

end
