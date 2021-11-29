% This script provides an inital guess for the design variables.
% The guess is data-informed (DI). We use experimental data to provide an
% initial guess of the joint variables but set constant values to the 
% muscle variable and the arm variables. We use a pre-defined final time 
% that is function of the imposed speed.
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function guess = getGuess_DI_tracking(Qs,nq,N,NMuscle,jointi,scaling)

%% Spline approximation of Qs to get Qdots and Qdotdots
Qs_spline.data = zeros(size(Qs.allinterpfilt));
Qs_spline.data(:,1) = Qs.allinterpfilt(:,1);
Qdots_spline.data = zeros(size(Qs.allinterpfilt));
Qdots_spline.data(:,1) = Qs.allinterpfilt(:,1);
Qdotdots_spline.data = zeros(size(Qs.allinterpfilt));
Qdotdots_spline.data(:,1) = Qs.allinterpfilt(:,1);
for i = 2:size(Qs.allinterpfilt,2)
    Qs.datafiltspline(i) = spline(Qs.allinterpfilt(:,1),Qs.allinterpfilt(:,i));
    [Qs_spline.data(:,i),Qdots_spline.data(:,i),...
        Qdotdots_spline.data(:,i)] = ...
        SplineEval_ppuval(Qs.datafiltspline(i),Qs.allinterpfilt(:,1),1);
end
% Filter the accelerations
order = 4;
cutoff_low = 10;
fs=1/mean(diff(Qs_spline.data(:,1)));
[af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
Qdotdots_spline.data(:,2:end) = filtfilt(af,bf,Qdotdots_spline.data(:,2:end));  

%% Qs: data-informed
% Pelvis tilt
guess.Qs(:,jointi.pelvis.tilt) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'));
% Pelvis list
guess.Qs(:,jointi.pelvis.list) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_list'));
% Pelvis rotation
guess.Qs(:,jointi.pelvis.rot) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rotation'));
% Pelvis_tx
guess.Qs(:,jointi.pelvis.tx) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
% Pelvis_ty
guess.Qs(:,jointi.pelvis.ty) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'));
% Pelvis_tz
guess.Qs(:,jointi.pelvis.tz) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz'));
% Hip flexion
guess.Qs(:,jointi.hip_flex.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'));
guess.Qs(:,jointi.hip_flex.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'));
% Hip adduction
guess.Qs(:,jointi.hip_add.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_l'));
guess.Qs(:,jointi.hip_add.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_r'));
% Hip rotation
guess.Qs(:,jointi.hip_rot.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_l'));
guess.Qs(:,jointi.hip_rot.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_r'));
% Knee
guess.Qs(:,jointi.knee.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'));
guess.Qs(:,jointi.knee.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'));
% Ankle
guess.Qs(:,jointi.ankle.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'));
guess.Qs(:,jointi.ankle.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'));
% Subtalar
guess.Qs(:,jointi.subt.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_l'));
guess.Qs(:,jointi.subt.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_r'));
% Trunk extension
guess.Qs(:,jointi.trunk.ext) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'));
% Trunk bending
guess.Qs(:,jointi.trunk.ben) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_bending'));
% Trunk rotation
guess.Qs(:,jointi.trunk.rot) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_rotation'));
% Arm flexion
guess.Qs(:,jointi.sh_flex.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_l'));
guess.Qs(:,jointi.sh_flex.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_r'));
% Arm adduction
guess.Qs(:,jointi.sh_add.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_l'));
guess.Qs(:,jointi.sh_add.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_r'));
% Arm rotation
guess.Qs(:,jointi.sh_rot.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'));
guess.Qs(:,jointi.sh_rot.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'));
% Elbow flexion
guess.Qs(:,jointi.elb.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'));
guess.Qs(:,jointi.elb.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'));

%% Qdots: data-informed
% Pelvis tilt
guess.Qdots(:,jointi.pelvis.tilt) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'));
% Pelvis list
guess.Qdots(:,jointi.pelvis.list) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_list'));
% Pelvis rotation
guess.Qdots(:,jointi.pelvis.rot) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rotation'));
% Pelvis_tx
guess.Qdots(:,jointi.pelvis.tx) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
% Pelvis_ty
guess.Qdots(:,jointi.pelvis.ty) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'));
% Pelvis_tz
guess.Qdots(:,jointi.pelvis.tz) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz'));
% Hip flexion
guess.Qdots(:,jointi.hip_flex.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'));
guess.Qdots(:,jointi.hip_flex.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'));
% Hip adduction
guess.Qdots(:,jointi.hip_add.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_l'));
guess.Qdots(:,jointi.hip_add.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_r'));
% Hip rotation
guess.Qdots(:,jointi.hip_rot.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_l'));
guess.Qdots(:,jointi.hip_rot.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_r'));
% Knee
guess.Qdots(:,jointi.knee.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'));
guess.Qdots(:,jointi.knee.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'));
% Ankle
guess.Qdots(:,jointi.ankle.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'));
guess.Qdots(:,jointi.ankle.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'));
% Subtalar
guess.Qdots(:,jointi.subt.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_l'));
guess.Qdots(:,jointi.subt.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_r'));
% Trunk extension
guess.Qdots(:,jointi.trunk.ext) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'));
% Trunk bending
guess.Qdots(:,jointi.trunk.ben) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_bending'));
% Trunk rotation
guess.Qdots(:,jointi.trunk.rot) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_rotation'));
% Arm flexion
guess.Qdots(:,jointi.sh_flex.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_l'));
guess.Qdots(:,jointi.sh_flex.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_r'));
% Arm adduction
guess.Qdots(:,jointi.sh_add.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_l'));
guess.Qdots(:,jointi.sh_add.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_r'));
% Arm rotation
guess.Qdots(:,jointi.sh_rot.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'));
guess.Qdots(:,jointi.sh_rot.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'));
% Elbow flexion
guess.Qdots(:,jointi.elb.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'));
guess.Qdots(:,jointi.elb.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'));

%% Qs and Qdots are intertwined
guess.QsQdots = zeros(N,2*nq.all);
guess.QsQdots(:,1:2:end) = guess.Qs;
guess.QsQdots(:,2:2:end) = guess.Qdots;

%% Qdotdots: data-informed
% Pelvis tilt
guess.Qdotdots(:,jointi.pelvis.tilt) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'));
% Pelvis list
guess.Qdotdots(:,jointi.pelvis.list) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_list'));
% Pelvis rotation
guess.Qdotdots(:,jointi.pelvis.rot) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_rotation'));
% Pelvis_tx
guess.Qdotdots(:,jointi.pelvis.tx) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
% Pelvis_ty
guess.Qdotdots(:,jointi.pelvis.ty) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'));
% Pelvis_tz
guess.Qdotdots(:,jointi.pelvis.tz) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tz'));
% Hip flexion
guess.Qdotdots(:,jointi.hip_flex.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'));
guess.Qdotdots(:,jointi.hip_flex.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'));
% Hip adduction
guess.Qdotdots(:,jointi.hip_add.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_l'));
guess.Qdotdots(:,jointi.hip_add.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_adduction_r'));
% Hip rotation
guess.Qdotdots(:,jointi.hip_rot.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_l'));
guess.Qdotdots(:,jointi.hip_rot.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rotation_r'));
% Knee
guess.Qdotdots(:,jointi.knee.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'));
guess.Qdotdots(:,jointi.knee.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'));
% Ankle
guess.Qdotdots(:,jointi.ankle.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'));
guess.Qdotdots(:,jointi.ankle.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'));
% Subtalar
guess.Qdotdots(:,jointi.subt.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_l'));
guess.Qdotdots(:,jointi.subt.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subtalar_angle_r'));
% Trunk extension
guess.Qdotdots(:,jointi.trunk.ext) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'));
% Trunk bending
guess.Qdotdots(:,jointi.trunk.ben) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_bending'));
% Trunk rotation
guess.Qdotdots(:,jointi.trunk.rot) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_rotation'));
% Arm flexion
guess.Qdotdots(:,jointi.sh_flex.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_l'));
guess.Qdotdots(:,jointi.sh_flex.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_flex_r'));
% Arm adduction
guess.Qdotdots(:,jointi.sh_add.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_l'));
guess.Qdotdots(:,jointi.sh_add.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_add_r'));
% Arm rotation
guess.Qdotdots(:,jointi.sh_rot.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_l'));
guess.Qdotdots(:,jointi.sh_rot.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'arm_rot_r'));
% Elbow flexion
guess.Qdotdots(:,jointi.elb.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_l'));
guess.Qdotdots(:,jointi.elb.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'elbow_flex_r'));

%% Muscle variables
guess.a = 0.1*ones(N,NMuscle);
guess.vA = 0.01*ones(N,NMuscle);
guess.FTtilde = 0.1*ones(N,NMuscle);
guess.dFTtilde = 0.01*ones(N,NMuscle);

%% Arm activations
guess.a_a = 0.1*ones(N,nq.arms);
guess.e_a = 0.1*ones(N,nq.arms);

%% Parameters contact model
% Original values
B_locSphere_s1_r    = [0.00190115788407966, -0.00382630379623308];
B_locSphere_s2_r    = [0.148386399942063, -0.028713422052654];
B_locSphere_s3_r    = [0.133001170607051, 0.0516362473449566];
B_locSphere_s4_r    = [0.06, -0.0187603084619177];    
B_locSphere_s5_r    = [0.0662346661991635, 0.0263641606741698];
B_locSphere_s6_r    = [0.045, 0.0618569567549652];
IG_rad              = 0.032*ones(1,6); 
guess.params = [B_locSphere_s1_r,B_locSphere_s2_r,...
    B_locSphere_s3_r,B_locSphere_s4_r,B_locSphere_s5_r,B_locSphere_s6_r,...
    IG_rad];

%% Scaling
guess.QsQdots   = guess.QsQdots./repmat(scaling.QsQdots,N,1);
guess.Qdotdots  = guess.Qdotdots./repmat(scaling.Qdotdots,N,1);
guess.a         = (guess.a)./repmat(scaling.a,N,size(guess.a,2));
guess.FTtilde   = (guess.FTtilde)./repmat(scaling.FTtilde,N,1);
guess.vA        = (guess.vA)./repmat(scaling.vA,N,size(guess.vA,2));
guess.dFTtilde  = (guess.dFTtilde)./repmat(scaling.dFTtilde,N,...
    size(guess.dFTtilde,2));
% no need to scale the IG for the arm activations / excitations
guess.params        = guess.params.*scaling.params.v + scaling.params.r;

end