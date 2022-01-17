% This script provides bounds and scaling factors for the design variables.
% The bounds on the joint variables are informed by experimental data.
% The bounds on the remaining variables are fixed.
% The bounds are scaled such that the upper/lower bounds cannot be
% larger/smaller than 1/-1.
%
% Author: Antoine Falisse
% Date: 12/19/2018
%
function [bounds,scaling] = getBounds_all_mtp(Qs,model_info,S)

nq = model_info.ExtFunIO.nq;
NMuscle = size(model_info.muscle_info.params.params,2);
joints = fields(model_info.ExtFunIO.coordi)';
jointi = model_info.ExtFunIO.coordi;

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
for i=[model_info.ExtFunIO.jointi.ground_pelvis model_info.ExtFunIO.jointi.back]
    bounds.Qs.upper(jointi.(joints{i})) = max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),joints{i}))));
    bounds.Qs.lower(jointi.(joints{i})) = min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),joints{i}))));
    bounds.Qdots.upper(jointi.(joints{i})) = max((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),joints{i}))));
    bounds.Qdots.lower(jointi.(joints{i})) = min((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),joints{i}))));
    bounds.Qdotdots.upper(jointi.(joints{i})) = max((Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),joints{i}))));
    bounds.Qdotdots.lower(jointi.(joints{i})) = min((Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),joints{i}))));
end
for i=[model_info.ExtFunIO.jointi.leg_r model_info.ExtFunIO.jointi.arm_r]
    bounds.Qs.upper(jointi.(joints{i})) = max(max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),joints{i}))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),[joints{i}(1:end-2) '_l']))));
    bounds.Qs.upper(jointi.([joints{i}(1:end-2) '_l'])) = bounds.Qs.upper(jointi.(joints{i}));
    bounds.Qs.lower(jointi.(joints{i})) = min(min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),joints{i}))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),[joints{i}(1:end-2) '_l']))));
    bounds.Qs.lower(jointi.([joints{i}(1:end-2) '_l'])) = bounds.Qs.lower(jointi.(joints{i}));

    bounds.Qdots.upper(jointi.(joints{i})) = max(max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),joints{i}))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),[joints{i}(1:end-2) '_l']))));
    bounds.Qdots.upper(jointi.([joints{i}(1:end-2) '_l'])) = bounds.Qdots.upper(jointi.(joints{i}));
    bounds.Qdots.lower(jointi.(joints{i})) = min(min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),joints{i}))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),[joints{i}(1:end-2) '_l']))));
    bounds.Qdots.lower(jointi.([joints{i}(1:end-2) '_l'])) = bounds.Qdots.lower(jointi.(joints{i}));

    bounds.Qdotdots.upper(jointi.(joints{i})) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),joints{i}))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),[joints{i}(1:end-2) '_l']))));
    bounds.Qdotdots.upper(jointi.([joints{i}(1:end-2) '_l'])) = bounds.Qdotdots.upper(jointi.(joints{i}));
    bounds.Qdotdots.lower(jointi.(joints{i})) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),joints{i}))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),[joints{i}(1:end-2) '_l']))));
    bounds.Qdotdots.lower(jointi.([joints{i}(1:end-2) '_l'])) = bounds.Qdotdots.lower(jointi.(joints{i}));
end
% The bounds are extended by twice the absolute difference between upper
% and lower bounds.
Qs_range = abs(bounds.Qs.upper - bounds.Qs.lower);
bounds.Qs.lower = bounds.Qs.lower - 2*Qs_range;
bounds.Qs.upper = bounds.Qs.upper + 2*Qs_range;
% For several joints, we manually adjust the bounds
% Pelvis_tx
bounds.Qs.upper(jointi.pelvis_tx) = 2;  
bounds.Qs.lower(jointi.pelvis_tx) = 0;
% Pelvis_ty
bounds.Qs.upper(jointi.pelvis_ty) = 1.1;  
bounds.Qs.lower(jointi.pelvis_ty) = 0.75;
% Pelvis_tz
bounds.Qs.upper(jointi.pelvis_tz) = 0.1;
bounds.Qs.lower(jointi.pelvis_tz) = -0.1;
% Mtp
bounds.Qs.upper(jointi.mtp_angle_l) = 1.05;
bounds.Qs.lower(jointi.mtp_angle_l) = -0.5;
bounds.Qs.upper(jointi.mtp_angle_r) = 1.05;
bounds.Qs.lower(jointi.mtp_angle_r) = -0.5;
% Elbow
bounds.Qs.lower(jointi.elbow_flex_l) = 0;
bounds.Qs.lower(jointi.elbow_flex_r) = 0;
% % Shoulder adduction
% bounds.Qs.upper(jointi.sh_add.l) = rec_uw_sh_add;
% bounds.Qs.upper(jointi.sh_add.r) = rec_uw_sh_add;
% % Shoulder rotation
% bounds.Qs.upper(jointi.sh_rot.l) = rec_uw_sh_rot;
% bounds.Qs.upper(jointi.sh_rot.r) = rec_uw_sh_rot;
% We adjust some bounds when we increase the speed to allow for the
% generation of running motions.
if S.v_tgt > 1.33
    % Pelvis tilt
    bounds.Qs.lower(jointi.pelvis_tilt) = -20*pi/180;
    % Shoulder flexion
    bounds.Qs.lower(jointi.arm_flex_l) = -50*pi/180;
    bounds.Qs.lower(jointi.arm_flex_r) = -50*pi/180;
end

% The bounds are extended by 3 times the absolute difference between upper
% and lower bounds.
Qdots_range = abs(bounds.Qdots.upper - bounds.Qdots.lower);
bounds.Qdots.lower = bounds.Qdots.lower - 3*Qdots_range;
bounds.Qdots.upper = bounds.Qdots.upper + 3*Qdots_range;
% Mtp
bounds.Qdots.upper(jointi.mtp_angle_l) = 13;
bounds.Qdots.lower(jointi.mtp_angle_l) = -13;
bounds.Qdots.upper(jointi.mtp_angle_r) = 13;
bounds.Qdots.lower(jointi.mtp_angle_r) = -13;
% We adjust some bounds when we increase the speed to allow for the
% generation of running motions.
if S.v_tgt > 1.33
    % Pelvis tx
    bounds.Qdots.upper(jointi.pelvis_tx) = 4;
end
% The bounds are extended by 3 times the absolute difference between upper
% and lower bounds.
Qdotdots_range = abs(bounds.Qdotdots.upper - bounds.Qdotdots.lower);
bounds.Qdotdots.lower = bounds.Qdotdots.lower - 3*Qdotdots_range;
bounds.Qdotdots.upper = bounds.Qdotdots.upper + 3*Qdotdots_range;
% Mtp
bounds.Qdotdots.upper(jointi.mtp_angle_l) = 500;
bounds.Qdotdots.lower(jointi.mtp_angle_l) = -500;
bounds.Qdotdots.upper(jointi.mtp_angle_r) = 500;
bounds.Qdotdots.lower(jointi.mtp_angle_r) = -500;

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

%% Mtp excitations
bounds.e_mtp.lower = -ones(1,nq.mtp);
bounds.e_mtp.upper = ones(1,nq.mtp);

%% Mtp activations
bounds.a_mtp.lower = -ones(1,nq.mtp);
bounds.a_mtp.upper = ones(1,nq.mtp);

%% Lumbar activations
% Only used when no muscles actuate the lumbar joints (e.g. Rajagopal
% model)
bounds.a_lumbar.lower = -ones(1,nq.trunk);
bounds.a_lumbar.upper = ones(1,nq.trunk);

%% Lumbar excitations
% Only used when no muscles actuate the lumbar joints (e.g. Rajagopal
% model)
bounds.e_lumbar.lower = -ones(1,nq.trunk);
bounds.e_lumbar.upper = ones(1,nq.trunk);

%% Final time
bounds.tf.lower = S.bounds.t_final.lower;
bounds.tf.upper = S.bounds.t_final.upper;

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
% Fixed scaling factor
scaling.LumbarTau = 150;
% Mtp torque actuators
% Fixed scaling factor
scaling.MtpTau = 100;
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
