% This script provides an inital guess for the design variables.
% The guess is quasi-random (QR). We set constant values to the muscle
% variables, the arm variables and most joint variables. We only ensure
% that the distance traveled is not null. The model is moving forward at a 
% constant speed and is standing on the ground. We use a pre-defined final
% time that is function of the imposed speed.
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function guess = getGuess_QR_opti(S,model_info,scaling,d)
N = S.solver.N_meshes; % number of mesh intervals
nq = model_info.ExtFunIO.jointi.nq;
NMuscle = model_info.muscle_info.NMuscle;
jointi = model_info.ExtFunIO.coordi;

%% Final time
% The final time is function of the imposed speed
all_speeds = 0.73:0.1:5;
all_tf = 0.70:-((0.70-0.35)/(length(all_speeds)-1)):0.35;
idx_speed = find(all_speeds==S.subject.v_pelvis_x_trgt);
if isempty(idx_speed)
    idx_speed = find(all_speeds > S.subject.v_pelvis_x_trgt,1,'first');
end
guess.tf = all_tf(idx_speed);

%% Qs
% The model is moving forward but with a standing position (Qs=0)
guess.Qs = zeros(N,nq.all);
guess.Qs(:,model_info.ExtFunIO.jointi.base_forward) = linspace(0,guess.tf*S.subject.v_pelvis_x_trgt,N);
% The model is standing on the ground
guess.Qs(:,jointi.pelvis_ty) = S.subject.IG_pelvis_y;

%% Qdots
guess.Qdots = zeros(N,nq.all);
% The model is moving forward with a constant speed
guess.Qdots(:,model_info.ExtFunIO.jointi.base_forward) = S.subject.v_pelvis_x_trgt;

%% Qdotdots
guess.Qdotdots = zeros(N,nq.all);

%% Muscle variables
guess.a = 0.1*ones(N,NMuscle);
guess.vA = 0.01*ones(N,NMuscle);
guess.FTtilde = 0.1*ones(N,NMuscle);
guess.dFTtilde = 0.01*ones(N,NMuscle);

%% Torque actuator activations
guess.a_a = 0.1*ones(N,nq.torqAct);
guess.e_a = 0.1*ones(N,nq.torqAct);

%% Add last mesh point to state variables
if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    % Lower limbs and trunk
    % Qs and Qdots are inverted after a half gait cycle BUT 
    % Pelvis: pelvis tilt and pelvis ty should be equal, pelvis
    % list, rot, tz should be opposite and pelvis tx should be equal
    % plus dist traveled.
    % Trunk: lumbar ext should be equal, lumbar bend and lumbar rot
    % should be of opposite.     
    % For "symmetric" joints, we invert right and left
    inv_X_Qs = zeros(1,nq.all);
    inv_X_Qs(1,model_info.ExtFunIO.symQs.QsInvA) = guess.Qs(1,model_info.ExtFunIO.symQs.QsInvB);
    inv_X_Qdots = zeros(1,nq.all);
    inv_X_Qdots(1,model_info.ExtFunIO.symQs.QdotsInvA) = guess.Qs(1,model_info.ExtFunIO.symQs.QdotsInvB);
    % For other joints, we take the opposite right and left
    inv_X_Qs(model_info.ExtFunIO.symQs.QsOpp) = ...
        -guess.Qs(1,model_info.ExtFunIO.symQs.QsOpp);           
    inv_X_Qdots(model_info.ExtFunIO.symQs.QsOpp) = ...
        -guess.Qdots(1,model_info.ExtFunIO.symQs.QsOpp);           
    dx = guess.Qs(end,jointi.pelvis_tx) - ...
        guess.Qs(end-1,jointi.pelvis_tx);
    inv_X_Qs(jointi.pelvis_tx) = ...
        guess.Qs(end,jointi.pelvis_tx) + dx;

    guess.Qs = [guess.Qs; inv_X_Qs];
    guess.Qdots = [guess.Qdots; inv_X_Qdots];
    guess.a = [guess.a; guess.a(1,model_info.ExtFunIO.symQs.MusInvB)];
    guess.FTtilde = [guess.FTtilde; guess.FTtilde(1,model_info.ExtFunIO.symQs.MusInvB)];
    guess.a_a = [guess.a_a; guess.a_a(1,:)];
else
    dx = guess.Qs(end,jointi.pelvis_tx) - ...
        guess.Qs(end-1,jointi.pelvis_tx);
    guess.Qs = [guess.Qs; guess.Qs(1,:)];
    guess.Qs(end,jointi.pelvis_tx) = guess.Qs(end-1,jointi.pelvis_tx) + dx;
    guess.Qdots = [guess.Qdots; guess.Qdots(1,:)];
    guess.a = [guess.a; guess.a(1,:)];
    guess.FTtilde = [guess.FTtilde; guess.FTtilde(1,:)];
    guess.a_a = [guess.a_a; guess.a_a(1,:)];
end


%% Scaling
guess.Qs = guess.Qs./repmat(scaling.Qs,N+1,1);
guess.Qdots = guess.Qdots./repmat(scaling.Qdots,N+1,1);
guess.Qdotdots  = guess.Qdotdots./repmat(scaling.Qdotdots,N,1);
guess.a         = (guess.a)./repmat(scaling.a,N+1,size(guess.a,2));
guess.FTtilde   = (guess.FTtilde)./repmat(scaling.FTtilde,N+1,1);
guess.vA        = (guess.vA)./repmat(scaling.vA,N,size(guess.vA,2));
guess.dFTtilde  = (guess.dFTtilde)./repmat(scaling.dFTtilde,N,size(guess.dFTtilde,2));
% guess.a_mtp_col = zeros(d*N,nq.mtp);
% guess.a_lumbar_col = zeros(d*N,nq.torso);

%% Collocation points
guess.a_col = zeros(d*N,NMuscle);
guess.FTtilde_col = zeros(d*N,NMuscle);
guess.Qs_col = zeros(d*N,nq.all);
guess.Qdots_col = zeros(d*N,nq.all);
guess.a_a_col = zeros(d*N,nq.arms);
guess.dFTtilde_col = zeros(d*N,NMuscle);
guess.Qdotdots_col = zeros(d*N,nq.all);
for k=1:N
    guess.a_col((k-1)*d+1:k*d,:) = repmat(guess.a(k,:),d,1); 
    guess.FTtilde_col((k-1)*d+1:k*d,:) = repmat(guess.FTtilde(k,:),d,1);
    guess.Qs_col((k-1)*d+1:k*d,:) = repmat(guess.Qs(k,:),d,1);
    guess.Qdots_col((k-1)*d+1:k*d,:) = repmat(guess.Qdots(k,:),d,1);
    guess.a_a_col((k-1)*d+1:k*d,:) = repmat(guess.a_a(k,:),d,1);
    guess.dFTtilde_col((k-1)*d+1:k*d,:) = repmat(guess.dFTtilde(k,:),d,1);
    guess.Qdotdots_col((k-1)*d+1:k*d,:) = repmat(guess.Qdotdots(k,:),d,1);
end
