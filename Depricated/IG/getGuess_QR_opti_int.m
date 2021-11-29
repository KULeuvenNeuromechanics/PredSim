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
function guess = getGuess_QR_opti_int(N,nq,NMuscle,scaling,v_tgt,jointi,d,PelvisY)


if ~exist('PelvisY','var')
    PelvisY = 0.9385;
end
%% Final time
% The final time is function of the imposed speed
all_speeds = 0.73:0.1:2.73;
all_tf = 0.70:-((0.70-0.35)/(length(all_speeds)-1)):0.35;
idx_speed = find(all_speeds==v_tgt);
if isempty(idx_speed)
    idx_speed = find(all_speeds > v_tgt,1,'first');
end
guess.tf = all_tf(idx_speed);

%% Qs
% The model is moving forward but with a standing position (Qs=0)
guess.Qs = zeros(N,nq.all);
guess.Qs(:,jointi.pelvis.tx) = linspace(0,guess.tf*v_tgt,N);
% The model is standing on the ground
guess.Qs(:,jointi.pelvis.ty) = PelvisY;

%% Qdots
guess.Qdots = zeros(N,nq.all);
% The model is moving forward with a constant speed
guess.Qdots(:,jointi.pelvis.tx) = v_tgt;
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
% Lower limbs and trunk
% Qs and Qdots are inverted after a half gait cycle BUT 
% Pelvis: pelvis tilt and pelvis ty should be equal, pelvis
% list, rot, tz should be opposite and pelvis tx should be equal
% plus dist traveled.
% Trunk: lumbar ext should be equal, lumbar bend and lumbar rot
% should be of opposite.     
orderQsInv = [jointi.pelvis.tilt:2*jointi.pelvis.tz,...
    2*jointi.hip_flex.r-1:2*jointi.hip_rot.r,...
    2*jointi.hip_flex.l-1:2*jointi.hip_rot.l,...
    2*jointi.knee.r-1:2*jointi.knee.r,...
    2*jointi.knee.l-1:2*jointi.knee.l,...
    2*jointi.ankle.r-1:2*jointi.ankle.r,...
    2*jointi.ankle.l-1:2*jointi.ankle.l,...
    2*jointi.subt.r-1:2*jointi.subt.r,...
    2*jointi.subt.l-1:2*jointi.subt.l,...
    2*jointi.mtp.r-1:2*jointi.mtp.r,...
    2*jointi.mtp.l-1:2*jointi.mtp.l,...
    2*jointi.trunk.ext-1:2*jointi.trunk.rot,...
    2*jointi.sh_flex.r-1:2*jointi.sh_rot.r,...
    2*jointi.sh_flex.l-1:2*jointi.sh_rot.l,...
    2*jointi.elb.r-1:2*jointi.elb.r,...
    2*jointi.elb.l-1:2*jointi.elb.l];        
orderQsOpp = [2*jointi.pelvis.list-1:2*jointi.pelvis.list,...   
    2*jointi.pelvis.rot-1:2*jointi.pelvis.rot,...
    2*jointi.pelvis.tz-1:2*jointi.pelvis.tz,...
    2*jointi.trunk.ben-1:2*jointi.trunk.ben,...
    2*jointi.trunk.rot-1:2*jointi.trunk.rot];
% For "symmetric" joints, we invert right and left
inv_X = guess.QsQdots(1,orderQsInv);
% For other joints, we take the opposite right and left
inv_X(orderQsOpp) = -guess.QsQdots(1,orderQsOpp);           
dx = guess.QsQdots(end,2*jointi.pelvis.tx-1) - ...
    guess.QsQdots(end-1,2*jointi.pelvis.tx-1);
inv_X(2*jointi.pelvis.tx-1) = ...
    guess.QsQdots(end,2*jointi.pelvis.tx-1) + dx;
guess.QsQdots = [guess.QsQdots; inv_X];

orderMusInv = [NMuscle/2+1:NMuscle,1:NMuscle/2];
guess.a = [guess.a; guess.a(1,orderMusInv)];
guess.FTtilde = [guess.FTtilde; guess.FTtilde(1,orderMusInv)];

orderArmInv = [jointi.sh_flex.r:jointi.sh_rot.r,...
    jointi.sh_flex.l:jointi.sh_rot.l,...
    jointi.elb.r:jointi.elb.r,...
    jointi.elb.l:jointi.elb.l]-jointi.sh_flex.l+1;
guess.a_a = [guess.a_a; guess.a_a(1,orderArmInv)];

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
