% This script provides bounds and scaling factors for the design variables.
% The bounds on the joint variables are informed by experimental data.
% The bounds on the remaining variables are fixed.
% The bounds are scaled such that the upper/lower bounds cannot be
% larger/smaller than 1/-1.
%
% Author: Antoine Falisse
% Date: 12/19/2018
%
function [bounds,scaling] = getBounds_tracking_mtp(Qs,NMuscle,nq,jointi,dev_cm,GRF,cs,csc)

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
% Filter the accelerations
order = 4;
cutoff_low = 10;
fs=1/mean(diff(Qs_spline.data(:,1)));
[af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
Qdotdots_spline.data(:,2:end) = filtfilt(af,bf,Qdotdots_spline.data(:,2:end));  

%% Qs: based on filtered experimental data (inverse kinematics)
% Pelvis tilt
bounds.Qs.upper(jointi.pelvis.tilt) = max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'))));
bounds.Qs.lower(jointi.pelvis.tilt) = min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'))));
% Pelvis list
bounds.Qs.upper(jointi.pelvis.list) = max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_list'))));
bounds.Qs.lower(jointi.pelvis.list) = min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_list'))));
% Pelvis rot
bounds.Qs.upper(jointi.pelvis.rot) = max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rotation'))));
bounds.Qs.lower(jointi.pelvis.rot) = min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rotation'))));
% Pelvis tx
bounds.Qs.upper(jointi.pelvis.tx) = max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'))));
bounds.Qs.lower(jointi.pelvis.tx) = min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'))));
% Pelvis ty
bounds.Qs.upper(jointi.pelvis.ty) = max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'))));
bounds.Qs.lower(jointi.pelvis.ty) = min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'))));
% Pelvis tz
bounds.Qs.upper(jointi.pelvis.tz) = max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz'))));
bounds.Qs.lower(jointi.pelvis.tz) = min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz'))));
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
% Mtp
bounds.Qs.upper(jointi.mtp.l) = max(max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'mtp_angle_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'mtp_angle_l'))));
bounds.Qs.lower(jointi.mtp.l) = min(min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'mtp_angle_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'mtp_angle_l'))));
bounds.Qs.upper(jointi.mtp.r) = bounds.Qs.upper(jointi.mtp.l);
bounds.Qs.lower(jointi.mtp.r) = bounds.Qs.lower(jointi.mtp.l);
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
% Shoulder rotation
bounds.Qs.upper(jointi.sh_rot.l) = max(max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'))));
bounds.Qs.lower(jointi.sh_rot.l) = min(min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'))));
bounds.Qs.upper(jointi.sh_rot.r) = bounds.Qs.upper(jointi.sh_rot.l);
bounds.Qs.lower(jointi.sh_rot.r) = bounds.Qs.lower(jointi.sh_rot.l);
% Elbow
bounds.Qs.upper(jointi.elb.l) = max(max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'))));
bounds.Qs.lower(jointi.elb.l) = min(min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'))));
bounds.Qs.upper(jointi.elb.r) = bounds.Qs.upper(jointi.elb.l);
bounds.Qs.lower(jointi.elb.r) = bounds.Qs.lower(jointi.elb.l);
% The bounds are extended by the absolute difference between upper
% and lower bounds.
Qs_range = abs(bounds.Qs.upper - bounds.Qs.lower);
bounds.Qs.lower = bounds.Qs.lower - Qs_range;
bounds.Qs.upper = bounds.Qs.upper + Qs_range;

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
% Mtp
bounds.Qdots.upper(jointi.mtp.l) = max(max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'mtp_angle_r'))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'mtp_angle_l'))));
bounds.Qdots.lower(jointi.mtp.l) = min(min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'mtp_angle_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'mtp_angle_l'))));
bounds.Qdots.upper(jointi.mtp.r) = bounds.Qdots.upper(jointi.mtp.l);
bounds.Qdots.lower(jointi.mtp.r) = bounds.Qdots.lower(jointi.mtp.l);
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
% The bounds are extended by the absolute difference between upper
% and lower bounds.
Qdots_range = abs(bounds.Qdots.upper - bounds.Qdots.lower);
bounds.Qdots.lower = bounds.Qdots.lower - Qdots_range;
bounds.Qdots.upper = bounds.Qdots.upper + Qdots_range;

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
% Mtp
bounds.Qdotdots.upper(jointi.mtp.l) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'mtp_angle_r'))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'mtp_angle_l'))));
bounds.Qdotdots.lower(jointi.mtp.l) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'mtp_angle_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'mtp_angle_l'))));
bounds.Qdotdots.upper(jointi.mtp.r) = bounds.Qdotdots.upper(jointi.mtp.l);
bounds.Qdotdots.lower(jointi.mtp.r) = bounds.Qdotdots.lower(jointi.mtp.l);
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
% The bounds are extended by the absolute difference between upper
% and lower bounds.
Qdotdots_range = abs(bounds.Qdotdots.upper - bounds.Qdotdots.lower);
bounds.Qdotdots.lower = bounds.Qdotdots.lower - Qdotdots_range;
bounds.Qdotdots.upper = bounds.Qdotdots.upper + Qdotdots_range;

%% Muscle activations
bounds.a.lower = zeros(1,NMuscle);
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

%% Mtp activations
bounds.a_mtp.lower = -ones(1,nq.mtp);
bounds.a_mtp.upper = ones(1,nq.mtp);

%% Mtp excitations
bounds.e_mtp.lower = -ones(1,nq.mtp);
bounds.e_mtp.upper = ones(1,nq.mtp);

%% Parameters contact model
% Original values & configuration (spheres 1,2,3&5 => calcn & 4&6 => toes)
B_locSphere_s1_r = [0.00190115788407966, -0.00382630379623308];
params_loc_IG = B_locSphere_s1_r;
B_locSphere_s2_r = [0.148386399942063, -0.028713422052654];
params_loc_IG = [params_loc_IG, B_locSphere_s2_r];
B_locSphere_s3_r = [0.133001170607051, 0.0516362473449566];
params_loc_IG = [params_loc_IG, B_locSphere_s3_r];
B_locSphere_s4_r = [0.06, -0.0187603084619177];    
params_loc_IG = [params_loc_IG, B_locSphere_s4_r];
if cs == 5 || cs == 6
    B_locSphere_s5_r = [0.0662346661991635, 0.0263641606741698];
    params_loc_IG = [params_loc_IG, B_locSphere_s5_r];
    if cs == 6
        B_locSphere_s6_r = [0.045, 0.0618569567549652];
        params_loc_IG = [params_loc_IG, B_locSphere_s6_r];
    end
end
radii = 0.032*ones(1,cs);
% Different configurations
if cs == 4 && csc == 2 
    % Sphere 3 => toes (originally sphere 4)
    % Sphere 4 => toes (originally sphere 6)
    % Original values
    params_loc_IG = B_locSphere_s1_r;
    params_loc_IG = [params_loc_IG, B_locSphere_s2_r];
    params_loc_IG = [params_loc_IG, B_locSphere_s4_r];
    B_locSphere_s6_r = [0.045, 0.0618569567549652];
    params_loc_IG = [params_loc_IG, B_locSphere_s6_r];
end
% Allowed deviations
params_loc_IG_upper = params_loc_IG + dev_cm.loc/1000;
params_loc_IG_lower = params_loc_IG - dev_cm.loc/1000;
radii_upper = radii + dev_cm.rad/100*radii;
radii_lower = radii - dev_cm.rad/100*radii;
% Manually constrained configurations
if cs == 4 && csc == 2
    params_loc_IG_lower(1) = 0;
    params_loc_IG_upper(3) = 0.15; 
    params_loc_IG_upper(4) = 0.03; 
    params_loc_IG_lower(4) = -0.03; 
    params_loc_IG_lower([5,7]) = 0.02;
    params_loc_IG_upper([6,8]) = 0.06;
    params_loc_IG_lower([6,8]) = -0.05;    
    radii_upper(1) = 0.065;
    % Case 12 vs case 11
    radii_upper(2) = 0.065;
end
if cs == 6 
    if csc == 2
        % Sphere 1
        B_locSphere_s1_r_dev_ub = B_locSphere_s1_r + 0.005;
        B_locSphere_s1_r_dev_lb = B_locSphere_s1_r - 0.005;
        radius_s1_ub = 0.032 + 0.05*0.032;
        radius_s1_lb = 0.032 - 0.05*0.032;
        % Sphere 5
        B_locSphere_s5_r_dev_ub = [B_locSphere_s5_r(1) + 0.01, B_locSphere_s5_r(2)];
        B_locSphere_s5_r_dev_lb = [B_locSphere_s5_r(1) - 0.01, B_locSphere_s5_r(2) - 0.01];
        radius_s5_ub = 0.032 + 0.30*0.032;
        radius_s5_lb = 0.032 - 0.30*0.032;
        % Sphere 3 
        B_locSphere_s3_r_dev_ub = [B_locSphere_s3_r(1), 0.04];
        B_locSphere_s3_r_dev_lb = [0.11, 0.02];
        radius_s3_ub = 0.032 + 0.30*0.032;
        radius_s3_lb = 0.032 - 0.30*0.032;
        % Sphere 2
        B_locSphere_s2_r_dev_ub = [0.155, -0.01];
        B_locSphere_s2_r_dev_lb = [0.11, B_locSphere_s2_r(2)];
        radius_s2_ub = 0.032 + 0.30*0.032;
        radius_s2_lb = 0.032 - 0.30*0.032;
        % Sphere 4
        B_locSphere_s4_r_dev_ub = [B_locSphere_s4_r(1), 0];
        B_locSphere_s4_r_dev_lb = [0, B_locSphere_s4_r(2)];
        radius_s4_ub = 0.032 + 0.30*0.032;
        radius_s4_lb = 0.032 - 0.30*0.032;
        % Sphere 6
        B_locSphere_s6_r_dev_ub = [0.025, 0.035];
        B_locSphere_s6_r_dev_lb = [0, 0];
        radius_s6_ub = 0.032 + 0.30*0.032;
        radius_s6_lb = 0.032 - 0.30*0.032;
    elseif csc == 3
        % Sphere 1
        B_locSphere_s1_r_dev_ub = B_locSphere_s1_r + 0.005;
        B_locSphere_s1_r_dev_lb = B_locSphere_s1_r - 0.005;
        radius_s1_ub = 0.032 + 0.05*0.032;
        radius_s1_lb = 0.032 - 0.05*0.032;
        % Sphere 5
        B_locSphere_s5_r_dev_ub = [B_locSphere_s5_r(1) + 0.01, B_locSphere_s5_r(2)];
        B_locSphere_s5_r_dev_lb = [B_locSphere_s5_r(1) - 0.01, B_locSphere_s5_r(2) - 0.01];
        radius_s5_ub = 0.032 + 0.05*0.032;
        radius_s5_lb = 0.032 - 0.05*0.032;
        % Sphere 3 
        B_locSphere_s3_r_dev_ub = [B_locSphere_s3_r(1), 0.04];
        B_locSphere_s3_r_dev_lb = [0.11, 0.02];
        radius_s3_ub = 0.032 + 0.50*0.032;
        radius_s3_lb = 0.032 - 0.50*0.032;
        % Sphere 2
        B_locSphere_s2_r_dev_ub = [0.155, -0.01];
        B_locSphere_s2_r_dev_lb = [0.11, B_locSphere_s2_r(2)];
        radius_s2_ub = 0.032 + 0.50*0.032;
        radius_s2_lb = 0.032 - 0.50*0.032;
        % Sphere 4
        B_locSphere_s4_r_dev_ub = [B_locSphere_s4_r(1), 0.035];
        B_locSphere_s4_r_dev_lb = [0, B_locSphere_s4_r(2)];
        radius_s4_ub = 0.032 + 0.50*0.032;
        radius_s4_lb = 0.032 - 0.50*0.032;
        % Sphere 6
        B_locSphere_s6_r_dev_ub = [B_locSphere_s4_r(1), 0.035];
        B_locSphere_s6_r_dev_lb = [0, B_locSphere_s4_r(2)];
        radius_s6_ub = 0.032 + 0.50*0.032;
        radius_s6_lb = 0.032 - 0.50*0.032;   
    elseif csc == 4
        % Sphere 1
        B_locSphere_s1_r_dev_ub = [0.005, 0.005];
        B_locSphere_s1_r_dev_lb = [-0.005, -0.005];
        radius_s1_ub = 0.032 + 0.01*0.032;
        radius_s1_lb = 0.032 - 0.50*0.032;
        % Sphere 2
        B_locSphere_s2_r_dev_ub = [0.07, 0.03];
        B_locSphere_s2_r_dev_lb = [-0.005, -0.02];
        radius_s2_ub = 0.032 + 0.01*0.032;
        radius_s2_lb = 0.032 - 0.50*0.032;
        % Sphere 3 
        B_locSphere_s3_r_dev_ub = [0.165, 0.04];
        B_locSphere_s3_r_dev_lb = [0.11, -0.03];
        radius_s3_ub = 0.032 + 0.50*0.032;
        radius_s3_lb = 0.032 - 0.50*0.032;
        % Sphere 4 
        B_locSphere_s4_r_dev_ub = [0.165, 0.04];
        B_locSphere_s4_r_dev_lb = [0.11, -0.03];
        radius_s4_ub = 0.032 + 0.50*0.032;
        radius_s4_lb = 0.032 - 0.50*0.032;
        % Sphere 5
        B_locSphere_s5_r_dev_ub = [0.06, 0.04];
        B_locSphere_s5_r_dev_lb = [0, -0.03];
        radius_s5_ub = 0.032 + 0.50*0.032;
        radius_s5_lb = 0.032 - 0.50*0.032;
        % Sphere 6
        B_locSphere_s6_r_dev_ub = [0.06, 0.04];
        B_locSphere_s6_r_dev_lb = [0, -0.03];
        radius_s6_ub = 0.032 + 0.50*0.032;
        radius_s6_lb = 0.032 - 0.50*0.032;
    elseif csc == 5
        % Sphere 1
        B_locSphere_s1_r_dev_ub = [0.005, 0.005];
        B_locSphere_s1_r_dev_lb = [-0.005, -0.005];
        radius_s1_ub = 0.032 + 0.01*0.032;
        radius_s1_lb = 0.032 - 0.50*0.032;
        % Sphere 2
        B_locSphere_s2_r_dev_ub = [0.07, 0.03];
        B_locSphere_s2_r_dev_lb = [-0.06, 0.02];
        radius_s2_ub = 0.032 + 0.01*0.032;
        radius_s2_lb = 0.032 - 0.50*0.032;
        % Sphere 3 
        B_locSphere_s3_r_dev_ub = [0.165, 0.04];
        B_locSphere_s3_r_dev_lb = [0.11, -0.03];
        radius_s3_ub = 0.032 + 0.50*0.032;
        radius_s3_lb = 0.032 - 0.50*0.032;
        % Sphere 4 
        B_locSphere_s4_r_dev_ub = [0.165, -0.01];
        B_locSphere_s4_r_dev_lb = [0.145, -0.03];
        radius_s4_ub = 0.032 + 0.50*0.032;
        radius_s4_lb = 0.032 - 0.50*0.032;
        % Sphere 5
        B_locSphere_s5_r_dev_ub = [0.06, 0.04];
        B_locSphere_s5_r_dev_lb = [0, -0.03];
        radius_s5_ub = 0.032 + 0.50*0.032;
        radius_s5_lb = 0.032 - 0.50*0.032;
        % Sphere 6
        B_locSphere_s6_r_dev_ub = [0.06, 0.04];
        B_locSphere_s6_r_dev_lb = [0, -0.03];
        radius_s6_ub = 0.032 + 0.50*0.032;
        radius_s6_lb = 0.032 - 0.50*0.032;
    elseif csc == 6
        % Sphere 1
        B_locSphere_s1_r_dev_ub = [0.005, 0.005];
        B_locSphere_s1_r_dev_lb = [-0.005, -0.005];
        radius_s1_ub = 0.032 + 0.01*0.032;
        radius_s1_lb = 0.032 - 0.90*0.032;
        % Sphere 2
        B_locSphere_s2_r_dev_ub = [0.07, 0.03];
        B_locSphere_s2_r_dev_lb = [0.06, 0.02];
        radius_s2_ub = 0.032 + 0.01*0.032;
        radius_s2_lb = 0.032 - 0.90*0.032;
        % Sphere 3 
        B_locSphere_s3_r_dev_ub = [0.165, 0.04];
        B_locSphere_s3_r_dev_lb = [0.11, -0.03];
        radius_s3_ub = 0.032 + 0.50*0.032;
        radius_s3_lb = 0.032 - 0.90*0.032;
        % Sphere 4 
        B_locSphere_s4_r_dev_ub = [0.165, -0.01];
        B_locSphere_s4_r_dev_lb = [0.145, -0.03];
        radius_s4_ub = 0.032 + 0.50*0.032;
        radius_s4_lb = 0.032 - 0.90*0.032;
        % Sphere 5
        B_locSphere_s5_r_dev_ub = [0.06, 0.04];
        B_locSphere_s5_r_dev_lb = [0, -0.03];
        radius_s5_ub = 0.032 + 0.50*0.032;
        radius_s5_lb = 0.032 - 0.90*0.032;
        % Sphere 6
        B_locSphere_s6_r_dev_ub = [0.06, 0.04];
        B_locSphere_s6_r_dev_lb = [0, -0.03];
        radius_s6_ub = 0.032 + 0.50*0.032;
        radius_s6_lb = 0.032 - 0.90*0.032;
    end    
    params_loc_IG_upper = [B_locSphere_s1_r_dev_ub,B_locSphere_s2_r_dev_ub,...
        B_locSphere_s3_r_dev_ub,B_locSphere_s4_r_dev_ub,...
        B_locSphere_s5_r_dev_ub,B_locSphere_s6_r_dev_ub];
    params_loc_IG_lower = [B_locSphere_s1_r_dev_lb,B_locSphere_s2_r_dev_lb,...
        B_locSphere_s3_r_dev_lb,B_locSphere_s4_r_dev_lb,...
        B_locSphere_s5_r_dev_lb,B_locSphere_s6_r_dev_lb];
    radii_upper = [radius_s1_ub,radius_s2_ub,radius_s3_ub,radius_s4_ub,...
        radius_s5_ub,radius_s6_ub];
    radii_lower = [radius_s1_lb,radius_s2_lb,radius_s3_lb,radius_s4_lb,...
        radius_s5_lb,radius_s6_lb];
end
if cs == 5
    if csc == 2
        % Sphere 1
        B_locSphere_s1_r_dev_ub = B_locSphere_s1_r + 0.005;
        B_locSphere_s1_r_dev_lb = B_locSphere_s1_r - 0.005;
        radius_s1_ub = 0.032 + 0.05*0.032;
        radius_s1_lb = 0.032 - 0.05*0.032;
        % Sphere 2
        B_locSphere_s2_r_dev_ub = [B_locSphere_s5_r(1) + 0.01, B_locSphere_s5_r(2)];
        B_locSphere_s2_r_dev_lb = [B_locSphere_s5_r(1) - 0.01, B_locSphere_s5_r(2) - 0.01];
        radius_s2_ub = 0.032 + 0.05*0.032;
        radius_s2_lb = 0.032 - 0.05*0.032;
        % Sphere 3 
        B_locSphere_s3_r_dev_ub = [B_locSphere_s3_r(1), 0.04];
        B_locSphere_s3_r_dev_lb = [0.11, 0.02];
        radius_s3_ub = 0.032 + 0.50*0.032;
        radius_s3_lb = 0.032 - 0.50*0.032;
        % Sphere 4
        B_locSphere_s4_r_dev_ub = [B_locSphere_s4_r(1), 0.035];
        B_locSphere_s4_r_dev_lb = [0, B_locSphere_s4_r(2)];
        radius_s4_ub = 0.032 + 0.50*0.032;
        radius_s4_lb = 0.032 - 0.50*0.032;
        % Sphere 5
        B_locSphere_s5_r_dev_ub = [B_locSphere_s4_r(1), 0.035];
        B_locSphere_s5_r_dev_lb = [0, B_locSphere_s4_r(2)];
        radius_s5_ub = 0.032 + 0.50*0.032;
        radius_s5_lb = 0.032 - 0.50*0.032;  
    elseif csc == 3
        % Sphere 1
        B_locSphere_s1_r_dev_ub = B_locSphere_s1_r + 0.005;
        B_locSphere_s1_r_dev_lb = B_locSphere_s1_r - 0.005;
        radius_s1_ub = 0.032 + 0.05*0.032;
        radius_s1_lb = 0.032 - 0.05*0.032;
        % Sphere 2
        B_locSphere_s2_r_dev_ub = [B_locSphere_s5_r(1) + 0.01, B_locSphere_s5_r(2)];
        B_locSphere_s2_r_dev_lb = [B_locSphere_s5_r(1) - 0.01, -0.02];
        radius_s2_ub = 0.032 + 0.05*0.032;
        radius_s2_lb = 0.032 - 0.05*0.032;
        % Sphere 3 
        B_locSphere_s3_r_dev_ub = [0.160, 0.04];
        B_locSphere_s3_r_dev_lb = [0.11, -0.03];
        radius_s3_ub = 0.032 + 0.50*0.032;
        radius_s3_lb = 0.032 - 0.50*0.032;
        % Sphere 4
        B_locSphere_s4_r_dev_ub = [B_locSphere_s4_r(1), 0.035];
        B_locSphere_s4_r_dev_lb = [0, B_locSphere_s4_r(2)];
        radius_s4_ub = 0.032 + 0.50*0.032;
        radius_s4_lb = 0.032 - 0.50*0.032;
        % Sphere 5
        B_locSphere_s5_r_dev_ub = [B_locSphere_s4_r(1), 0.035];
        B_locSphere_s5_r_dev_lb = [0, B_locSphere_s4_r(2)];
        radius_s5_ub = 0.032 + 0.50*0.032;
        radius_s5_lb = 0.032 - 0.50*0.032; 
    end    
    params_loc_IG_upper = [B_locSphere_s1_r_dev_ub,B_locSphere_s2_r_dev_ub,...
        B_locSphere_s3_r_dev_ub,B_locSphere_s4_r_dev_ub,...
        B_locSphere_s5_r_dev_ub];
    params_loc_IG_lower = [B_locSphere_s1_r_dev_lb,B_locSphere_s2_r_dev_lb,...
        B_locSphere_s3_r_dev_lb,B_locSphere_s4_r_dev_lb,...
        B_locSphere_s5_r_dev_lb];
    radii_upper = [radius_s1_ub,radius_s2_ub,radius_s3_ub,radius_s4_ub,...
        radius_s5_ub];
    radii_lower = [radius_s1_lb,radius_s2_lb,radius_s3_lb,radius_s4_lb,...
        radius_s5_lb];
end
if cs == 4 && csc == 3
    % Sphere 1
    B_locSphere_s1_r_dev_ub = B_locSphere_s1_r + 0.005;
    B_locSphere_s1_r_dev_lb = B_locSphere_s1_r - 0.005;
    radius_s1_ub = 0.032 + 0.05*0.032;
    radius_s1_lb = 0.032 - 0.05*0.032;
    % Sphere 2 
    B_locSphere_s2_r_dev_ub = [0.160, 0.04];
    B_locSphere_s2_r_dev_lb = [0.11, -0.03];
    radius_s2_ub = 0.032 + 0.50*0.032;
    radius_s2_lb = 0.032 - 0.50*0.032;
    % Sphere 3
    B_locSphere_s3_r_dev_ub = [B_locSphere_s4_r(1), 0.035];
    B_locSphere_s3_r_dev_lb = [0, B_locSphere_s4_r(2)];
    radius_s3_ub = 0.032 + 0.50*0.032;
    radius_s3_lb = 0.032 - 0.50*0.032;
    % Sphere 4
    B_locSphere_s4_r_dev_ub = [B_locSphere_s4_r(1), 0.035];
    B_locSphere_s4_r_dev_lb = [0, B_locSphere_s4_r(2)];
    radius_s4_ub = 0.032 + 0.50*0.032;
    radius_s4_lb = 0.032 - 0.50*0.032; 
    params_loc_IG_upper = [B_locSphere_s1_r_dev_ub,B_locSphere_s2_r_dev_ub,...
        B_locSphere_s3_r_dev_ub,B_locSphere_s4_r_dev_ub];
    params_loc_IG_lower = [B_locSphere_s1_r_dev_lb,B_locSphere_s2_r_dev_lb,...
        B_locSphere_s3_r_dev_lb,B_locSphere_s4_r_dev_lb];
    radii_upper = [radius_s1_ub,radius_s2_ub,radius_s3_ub,radius_s4_ub];
    radii_lower = [radius_s1_lb,radius_s2_lb,radius_s3_lb,radius_s4_lb];
end
bounds.params.upper = [params_loc_IG_upper,radii_upper];
bounds.params.lower = [params_loc_IG_lower,radii_lower];

%% GRF
bounds.GRF.lower = min(GRF.val.all(:,2:end));
bounds.GRF.upper = max(GRF.val.all(:,2:end));
% Extend bounds to give some flexibility
GRF_ROM = abs(bounds.GRF.upper - bounds.GRF.lower);
bounds.GRF.lower = bounds.GRF.lower - GRF_ROM;
bounds.GRF.upper = bounds.GRF.upper + GRF_ROM;

%% GRM
bounds.GRM.lower = min(GRF.MorGF.allinterp(:,2:end));
bounds.GRM.upper = max(GRF.MorGF.allinterp(:,2:end));
% Extend bounds to give some flexibility
GRM_ROM = abs(bounds.GRM.upper - bounds.GRM.lower);
bounds.GRM.lower = bounds.GRM.lower - GRM_ROM;
bounds.GRM.upper = bounds.GRM.upper + GRM_ROM;

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
% Joint torques
% Fixed scaling factor
Topt = 150;
scaling.T = repmat(Topt,1,nq.act);
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
% Mtp activations
scaling.a_mtp = 1;
% Mtp excitations
scaling.e_mtp = 1;
% Fixed scaling factor
scaling.MtpTau = 100;
% Time derivative of muscle-tendon forces
% Fixed scaling factor
scaling.dFTtilde = 100;
% Muscle-tendon forces
scaling.FTtilde         = max(...
    abs(bounds.FTtilde.lower),abs(bounds.FTtilde.upper)); 
bounds.FTtilde.lower    = (bounds.FTtilde.lower)./scaling.FTtilde;
bounds.FTtilde.upper    = (bounds.FTtilde.upper)./scaling.FTtilde;
% Parameters contact model
scaling.params.v = 1./(bounds.params.upper-bounds.params.lower);
scaling.params.r = ...
    0.5 - bounds.params.upper./(bounds.params.upper-bounds.params.lower);
bounds.params.lower = -0.5*ones(1,length(scaling.params.v));
bounds.params.upper = 0.5*ones(1,length(scaling.params.v));
% GRF
scaling.GRF = max(abs(bounds.GRF.lower),abs(bounds.GRF.upper));
bounds.GRF.lower = (bounds.GRF.lower)./scaling.GRF;
bounds.GRF.upper = (bounds.GRF.upper)./scaling.GRF;
% GRM
scaling.GRM      = max(abs(bounds.GRM.lower),abs(bounds.GRM.upper));
bounds.GRM.lower = (bounds.GRM.lower)./scaling.GRM;
bounds.GRM.upper = (bounds.GRM.upper)./scaling.GRM;

end
