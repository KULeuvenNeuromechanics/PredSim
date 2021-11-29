% This script provides an inital guess for the design variables.
% The guess is data-informed (DI). We use experimental data to provide an
% initial guess of the joint variables but set constant values to the 
% muscle variable and the arm variables. We use a pre-defined final time 
% that is function of the imposed speed.
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function guess = getGuess_DI(Qs,nq,N,time_IC,NMuscle,jointi,scaling,v_tgt)

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

%% Qs: data-informed
% Pelvis tilt
guess.Qs_all.data(:,jointi.pelvis.tilt) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'));
% Pelvis list
guess.Qs_all.data(:,jointi.pelvis.list) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_list'));
% Pelvis rotation
guess.Qs_all.data(:,jointi.pelvis.rot) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rotation'));
% Pelvis_tx
guess.Qs_all.data(:,jointi.pelvis.tx) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
% Pelvis_ty
guess.Qs_all.data(:,jointi.pelvis.ty) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'));
% Pelvis_tz
guess.Qs_all.data(:,jointi.pelvis.tz) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz'));
% Hip flexion
guess.Qs_all.data(:,jointi.hip_flex.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'));
guess.Qs_all.data(:,jointi.hip_flex.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'));
% Hip adduction
guess.Qs_all.data(:,jointi.hip_add.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_l'));
guess.Qs_all.data(:,jointi.hip_add.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_r'));
% Hip rotation
guess.Qs_all.data(:,jointi.hip_rot.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_l'));
guess.Qs_all.data(:,jointi.hip_rot.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_r'));
% Knee
guess.Qs_all.data(:,jointi.knee.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'));
guess.Qs_all.data(:,jointi.knee.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'));
% Ankle
guess.Qs_all.data(:,jointi.ankle.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'));
guess.Qs_all.data(:,jointi.ankle.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'));
% Subtalar
guess.Qs_all.data(:,jointi.subt.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_l'));
guess.Qs_all.data(:,jointi.subt.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_r'));
% Trunk extension
guess.Qs_all.data(:,jointi.trunk.ext) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'));
% Trunk bending
guess.Qs_all.data(:,jointi.trunk.ben) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_bending'));
% Trunk rotation
guess.Qs_all.data(:,jointi.trunk.rot) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_rotation'));
% Arm flexion
guess.Qs_all.data(:,jointi.sh_flex.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_l'));
guess.Qs_all.data(:,jointi.sh_flex.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_r'));
% Arm adduction
guess.Qs_all.data(:,jointi.sh_add.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_l'));
guess.Qs_all.data(:,jointi.sh_add.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_r'));
% Arm rotation
guess.Qs_all.data(:,jointi.sh_rot.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'));
guess.Qs_all.data(:,jointi.sh_rot.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'));
% Elbow flexion
guess.Qs_all.data(:,jointi.elb.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'));
guess.Qs_all.data(:,jointi.elb.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'));
% Interpolation
Qs_time = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'time'));
time_expi.Qs(1) = find(round(Qs_time,3) == round(time_IC(1),3));
time_expi.Qs(2) = find(round(Qs_time,3) == round(time_IC(2),3));
step = (Qs_time(time_expi.Qs(2))-Qs_time(time_expi.Qs(1)))/(N-1);
interval = Qs_time(time_expi.Qs(1)):step:Qs_time(time_expi.Qs(2));
guess.Qs = interp1(round(Qs_time,4),guess.Qs_all.data,round(interval,4));
guess.Qs(:,jointi.pelvis.tx) = guess.Qs(:,jointi.pelvis.tx) - ....
    guess.Qs(1,jointi.pelvis.tx);

%% Qdots: data-informed
% Pelvis tilt
guess.Qdots_all.data(:,jointi.pelvis.tilt) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'));
% Pelvis list
guess.Qdots_all.data(:,jointi.pelvis.list) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_list'));
% Pelvis rotation
guess.Qdots_all.data(:,jointi.pelvis.rot) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rotation'));
% Pelvis_tx
guess.Qdots_all.data(:,jointi.pelvis.tx) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
% Pelvis_ty
guess.Qdots_all.data(:,jointi.pelvis.ty) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'));
% Pelvis_tz
guess.Qdots_all.data(:,jointi.pelvis.tz) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz'));
% Hip flexion
guess.Qdots_all.data(:,jointi.hip_flex.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'));
guess.Qdots_all.data(:,jointi.hip_flex.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'));
% Hip adduction
guess.Qdots_all.data(:,jointi.hip_add.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_l'));
guess.Qdots_all.data(:,jointi.hip_add.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_r'));
% Hip rotation
guess.Qdots_all.data(:,jointi.hip_rot.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_l'));
guess.Qdots_all.data(:,jointi.hip_rot.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_r'));
% Knee
guess.Qdots_all.data(:,jointi.knee.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'));
guess.Qdots_all.data(:,jointi.knee.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'));
% Ankle
guess.Qdots_all.data(:,jointi.ankle.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'));
guess.Qdots_all.data(:,jointi.ankle.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'));
% Subtalar
guess.Qdots_all.data(:,jointi.subt.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_l'));
guess.Qdots_all.data(:,jointi.subt.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_r'));
% Trunk extension
guess.Qdots_all.data(:,jointi.trunk.ext) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'));
% Trunk bending
guess.Qdots_all.data(:,jointi.trunk.ben) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_bending'));
% Trunk rotation
guess.Qdots_all.data(:,jointi.trunk.rot) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_rotation'));
% Arm flexion
guess.Qdots_all.data(:,jointi.sh_flex.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_l'));
guess.Qdots_all.data(:,jointi.sh_flex.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_r'));
% Arm adduction
guess.Qdots_all.data(:,jointi.sh_add.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_l'));
guess.Qdots_all.data(:,jointi.sh_add.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_r'));
% Arm rotation
guess.Qdots_all.data(:,jointi.sh_rot.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'));
guess.Qdots_all.data(:,jointi.sh_rot.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'));
% Elbow flexion
guess.Qdots_all.data(:,jointi.elb.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'));
guess.Qdots_all.data(:,jointi.elb.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'));
% Interpolation
guess.Qdots = interp1(round(Qs_time,4),guess.Qdots_all.data,...
    round(interval,4));

% Qs and Qdots are intertwined
guess.QsQdots = zeros(N,2*nq.all);
guess.QsQdots(:,1:2:end) = guess.Qs;
guess.QsQdots(:,2:2:end) = guess.Qdots;

%% Qdotdots: data-informed
% Pelvis tilt
guess.Qdotdots_all.data(:,jointi.pelvis.tilt) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'));
% Pelvis list
guess.Qdotdots_all.data(:,jointi.pelvis.list) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_list'));
% Pelvis rotation
guess.Qdotdots_all.data(:,jointi.pelvis.rot) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rotation'));
% Pelvis_tx
guess.Qdotdots_all.data(:,jointi.pelvis.tx) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
% Pelvis_ty
guess.Qdotdots_all.data(:,jointi.pelvis.ty) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'));
% Pelvis_tz
guess.Qdotdots_all.data(:,jointi.pelvis.tz) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz'));
% Hip flexion
guess.Qdotdots_all.data(:,jointi.hip_flex.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'));
guess.Qdotdots_all.data(:,jointi.hip_flex.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'));
% Hip adduction
guess.Qdotdots_all.data(:,jointi.hip_add.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_l'));
guess.Qdotdots_all.data(:,jointi.hip_add.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_r'));
% Hip rotation
guess.Qdotdots_all.data(:,jointi.hip_rot.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_l'));
guess.Qdotdots_all.data(:,jointi.hip_rot.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_r'));
% Knee
guess.Qdotdots_all.data(:,jointi.knee.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'));
guess.Qdotdots_all.data(:,jointi.knee.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'));
% Ankle
guess.Qdotdots_all.data(:,jointi.ankle.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'));
guess.Qdotdots_all.data(:,jointi.ankle.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'));
% Subtalar
guess.Qdotdots_all.data(:,jointi.subt.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_l'));
guess.Qdotdots_all.data(:,jointi.subt.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_r'));
% Trunk extension
guess.Qdotdots_all.data(:,jointi.trunk.ext) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'));
% Trunk bending
guess.Qdotdots_all.data(:,jointi.trunk.ben) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_bending'));
% Trunk rotation
guess.Qdotdots_all.data(:,jointi.trunk.rot) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_rotation'));
% Arm flexion
guess.Qdotdots_all.data(:,jointi.sh_flex.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_l'));
guess.Qdotdots_all.data(:,jointi.sh_flex.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_r'));
% Arm adduction
guess.Qdotdots_all.data(:,jointi.sh_add.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_l'));
guess.Qdotdots_all.data(:,jointi.sh_add.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_r'));
% Arm rotation
guess.Qdotdots_all.data(:,jointi.sh_rot.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'));
guess.Qdotdots_all.data(:,jointi.sh_rot.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'));
% Elbow flexion
guess.Qdotdots_all.data(:,jointi.elb.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'));
guess.Qdotdots_all.data(:,jointi.elb.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'));
% Interpolation
guess.Qdotdots = interp1(round(Qs_time,4),guess.Qdotdots_all.data,...
    round(interval,4));

%% Muscle variables
guess.a = 0.1*ones(N,NMuscle);
guess.vA = 0.01*ones(N,NMuscle);
guess.FTtilde = 0.1*ones(N,NMuscle);
guess.dFTtilde = 0.01*ones(N,NMuscle);

%% Arm activations
guess.a_a = 0.1*ones(N,nq.arms);
guess.e_a = 0.1*ones(N,nq.arms);

%% Final time
% The final time is function of the imposed speed
all_speeds = 0.73:0.1:2.73;
all_tf = 0.70:-((0.70-0.35)/(length(all_speeds)-1)):0.35;
idx_speed = find(all_speeds==v_tgt);
if isempty(idx_speed)
    idx_speed = find(all_speeds > v_tgt,1,'first');
end
guess.tf = all_tf(idx_speed);

%% Scaling
guess.QsQdots = guess.QsQdots./repmat(scaling.QsQdots,N,1);
guess.Qdotdots = guess.Qdotdots./repmat(scaling.Qdotdots,N,1);
guess.a         = (guess.a)./repmat(scaling.a,N,size(guess.a,2));
guess.FTtilde   = (guess.FTtilde)./repmat(scaling.FTtilde,N,1);
guess.vA        = (guess.vA)./repmat(scaling.vA,N,size(guess.vA,2));
guess.dFTtilde  = (guess.dFTtilde)./repmat(scaling.dFTtilde,N,...
    size(guess.dFTtilde,2));

end
