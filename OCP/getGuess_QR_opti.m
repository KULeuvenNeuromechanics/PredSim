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
function guess = getGuess_QR_opti(N,scaling,model_info,S,d)
nq = model_info.ExtFunIO.nq;
% NMuscle = size(model_info.muscle_info.params.params,2);
NMuscle = length(model_info.muscle_info.muscle_names);
% joints = fields(model_info.ExtFunIO.coordi)';
jointi = model_info.ExtFunIO.coordi;
PelvisY = S.subject.IG_pelvis_y;

if ~exist('PelvisY','var')
    PelvisY = 0.9385;
end
%% Final time
% The final time is function of the imposed speed
all_speeds = 0.73:0.1:2.73;
all_tf = 0.70:-((0.70-0.35)/(length(all_speeds)-1)):0.35;
idx_speed = find(all_speeds==S.subject.v_pelvis_x_trgt);
if isempty(idx_speed)
    idx_speed = find(all_speeds > S.subject.v_pelvis_x_trgt,1,'first');
end
guess.tf = all_tf(idx_speed);

%% Qs
% The model is moving forward but with a standing position (Qs=0)
guess.Qs = zeros(N,nq.all);
guess.Qs(:,jointi.pelvis_tx) = linspace(0,guess.tf*S.subject.v_pelvis_x_trgt,N);
% The model is standing on the ground
guess.Qs(:,jointi.pelvis_ty) = PelvisY;

%% Qdots
guess.Qdots = zeros(N,nq.all);
% The model is moving forward with a constant speed
guess.Qdots(:,jointi.pelvis_tx) = S.subject.v_pelvis_x_trgt;
% Qs and Qdots are intertwined
guess.QsQdots = zeros(N,2*nq.all);
guess.QsQdots(:,1:2:end) = guess.Qs;
guess.QsQdots(:,2:2:end) = guess.Qdots;

%% Qdotdots
guess.Qdotdots = zeros(N,nq.all);

%% Muscle variables
guess.a = 0.1*ones(N,NMuscle);
guess.vA = 0.01*ones(N,NMuscle);
guess.FTtilde = 0.1*ones(N,NMuscle);
guess.dFTtilde = 0.01*ones(N,NMuscle);

%% Arm activations
guess.a_a = 0.1*ones(N,nq.arms);
guess.e_a = 0.1*ones(N,nq.arms);

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
    inv_X = guess.QsQdots(1,model_info.ExtFunIO.symQs.orderQsInv);
    % For other joints, we take the opposite right and left
    inv_X(model_info.ExtFunIO.symQs.orderQsOpp1) = ...
        -guess.QsQdots(1,model_info.ExtFunIO.symQs.orderQsOpp1);           
    dx = guess.QsQdots(end,2*jointi.pelvis_tx-1) - ...
        guess.QsQdots(end-1,2*jointi.pelvis_tx-1);
    inv_X(2*jointi.pelvis.tx-1) = ...
        guess.QsQdots(end,2*jointi.pelvis_tx-1) + dx;

    orderMusInv = [NMuscle/2+1:NMuscle,1:NMuscle/2];

    guess.QsQdots = [guess.QsQdots; inv_X];
    guess.a = [guess.a; guess.a(1,orderMusInv)];
    guess.FTtilde = [guess.FTtilde; guess.FTtilde(1,orderMusInv)];
    guess.a_a = [guess.a_a; guess.a_a(1,model_info.ExtFunIO.symQs.orderArmInv)];
else
    guess.QsQdots = [guess.QsQdots; guess.QsQdots(1,:)];
    guess.a = [guess.a; guess.a(1,:)];
    guess.FTtilde = [guess.FTtilde; guess.FTtilde(1,:)];
    guess.a_a = [guess.a_a; guess.a_a(1,:)];
end

%% Mtp activations
guess.a_mtp = 0.1*ones(N+1,nq.mtp);
guess.e_mtp = 0.1*ones(N,nq.mtp);

%% Mtp lumbar activations
% Only used when no muscles actuate the lumbar joints (e.g. Rajagopal
% model)
guess.a_lumbar = 0.1*ones(N+1,nq.trunk);
guess.e_lumbar = 0.1*ones(N,nq.trunk);

%% Scaling
guess.QsQdots   = guess.QsQdots./repmat(scaling.QsQdots,N+1,1);
guess.Qdotdots  = guess.Qdotdots./repmat(scaling.Qdotdots,N,1);
guess.a         = (guess.a)./repmat(scaling.a,N+1,size(guess.a,2));
guess.FTtilde   = (guess.FTtilde)./repmat(scaling.FTtilde,N+1,1);
guess.vA        = (guess.vA)./repmat(scaling.vA,N,size(guess.vA,2));
guess.dFTtilde  = (guess.dFTtilde)./repmat(scaling.dFTtilde,N,...
    size(guess.dFTtilde,2));
guess.a_mtp_col = zeros(d*N,nq.mtp);
guess.a_lumbar_col = zeros(d*N,nq.trunk);

%% Collocation points
    guess.a_col = zeros(d*N,NMuscle);
    guess.FTtilde_col = zeros(d*N,NMuscle);
    guess.QsQdots_col = zeros(d*N,2*nq.all);
    guess.a_a_col = zeros(d*N,nq.arms);
    guess.dFTtilde_col = zeros(d*N,NMuscle);
    guess.Qdotdots_col = zeros(d*N,nq.all);
for k=1:N
    guess.a_col((k-1)*d+1:k*d,:) = repmat(guess.a(k,:),d,1); 
    guess.FTtilde_col((k-1)*d+1:k*d,:) = repmat(guess.FTtilde(k,:),d,1);
    guess.QsQdots_col((k-1)*d+1:k*d,:) = repmat(guess.QsQdots(k,:),d,1);
    guess.a_a_col((k-1)*d+1:k*d,:) = repmat(guess.a_a(k,:),d,1);
    guess.dFTtilde_col((k-1)*d+1:k*d,:) = repmat(guess.dFTtilde(k,:),d,1);
    guess.Qdotdots_col((k-1)*d+1:k*d,:) = repmat(guess.Qdotdots(k,:),d,1);
end
