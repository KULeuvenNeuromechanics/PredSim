% This script provides an inital guess for the design variables.
% The guess is data-informed (DI). We use experimental data to provide an
% initial guess of the joint variables but set constant values to the 
% muscle variable and the arm variables. We use a pre-defined final time 
% that is function of the imposed speed.
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function guess = getGuess_DI_opti(Qs,N,time_IC,scaling,S,d,model_info)
nq = model_info.ExtFunIO.nq;
% NMuscle = size(model_info.muscle_info.params.params,2);
NMuscle = length(model_info.muscle_info.muscle_names);
coordinate_names = fieldnames(model_info.ExtFunIO.coordi);
NCoord = length(coordinate_names);
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

for i=1:NCoord
    coordinate = coordinate_names{i};
    coord_idx = model_info.ExtFunIO.coordi.(coordinate);
    guess.Qs_all.data(:,coord_idx) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),coordinate));
    guess.Qdots_all.data(:,coord_idx) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),coordinate));
    guess.Qdotdots_all.data(:,coord_idx) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),coordinate));
    if contains(coordinate,'mtp')
        guess.Qs_all.data(:,coord_idx) = zeros(size(guess.Qs_all.data,1),1);
        guess.Qdots_all.data(:,coord_idx) = zeros(size(guess.Qdots_all.data,1),1);
        guess.Qdotdots_all.data(:,coord_idx) = zeros(size(guess.Qdotdots_all.data,1),1);
    end
end
    
% Interpolation
Qs_time = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'time'));
time_expi.Qs(1) = find(round(Qs_time,3) == round(time_IC(1),3));
time_expi.Qs(2) = find(round(Qs_time,3) == round(time_IC(2),3));
step = (Qs_time(time_expi.Qs(2))-Qs_time(time_expi.Qs(1)))/(N-1);
interval = Qs_time(time_expi.Qs(1)):step:Qs_time(time_expi.Qs(2));
guess.Qs = interp1(round(Qs_time,4),guess.Qs_all.data,round(interval,4));
guess.Qs(:,jointi.pelvis_tx) = guess.Qs(:,jointi.pelvis_tx) - ....
    guess.Qs(1,jointi.pelvis_tx);

% Interpolation
guess.Qdots = interp1(round(Qs_time,4),guess.Qdots_all.data,...
    round(interval,4));

% Qs and Qdots are intertwined
guess.QsQdots = zeros(N,2*nq.all);
guess.QsQdots(:,1:2:end) = guess.Qs;
guess.QsQdots(:,2:2:end) = guess.Qdots;

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

%% Final time
% The final time is function of the imposed speed
all_speeds = 0.73:0.1:5;
all_tf = 0.70:-((0.70-0.35)/(length(all_speeds)-1)):0.35;
idx_speed = find(all_speeds==S.subject.vPelvis_x_trgt);
if isempty(idx_speed)
    idx_speed = find(all_speeds > S.subject.vPelvis_x_trgt,1,'first');
end
guess.tf = all_tf(idx_speed);

%% Scaling
guess.QsQdots = guess.QsQdots./repmat(scaling.QsQdots,N+1,1);
guess.Qdotdots = guess.Qdotdots./repmat(scaling.Qdotdots,N,1);
guess.a         = (guess.a)./repmat(scaling.a,N+1,size(guess.a,2));
guess.FTtilde   = (guess.FTtilde)./repmat(scaling.FTtilde,N+1,1);
guess.vA        = (guess.vA)./repmat(scaling.vA,N,size(guess.vA,2));
guess.dFTtilde  = (guess.dFTtilde)./repmat(scaling.dFTtilde,N,...
    size(guess.dFTtilde,2));


%% Collocation points
guess.a_col = zeros(d*N,NMuscle);
guess.FTtilde_col = zeros(d*N,NMuscle);
guess.QsQdots_col = zeros(d*N,2*nq.all);
guess.a_a_col = zeros(d*N,nq.arms);
guess.a_mtp_col = zeros(d*N,nq.mtp);
guess.dFTtilde_col = zeros(d*N,NMuscle);
guess.Qdotdots_col = zeros(d*N,nq.all);
guess.a_lumbar_col = zeros(d*N,nq.trunk);
for k=1:N
    guess.a_col((k-1)*d+1:k*d,:) = repmat(guess.a(k,:),d,1); 
    guess.FTtilde_col((k-1)*d+1:k*d,:) = repmat(guess.FTtilde(k,:),d,1);
    guess.QsQdots_col((k-1)*d+1:k*d,:) = repmat(guess.QsQdots(k,:),d,1);
    guess.a_a_col((k-1)*d+1:k*d,:) = repmat(guess.a_a(k,:),d,1);
    guess.a_mtp_col((k-1)*d+1:k*d,:) = repmat(guess.a_mtp(k,:),d,1);
    guess.dFTtilde_col((k-1)*d+1:k*d,:) = repmat(guess.dFTtilde(k,:),d,1);
    guess.Qdotdots_col((k-1)*d+1:k*d,:) = repmat(guess.Qdotdots(k,:),d,1);
end
end
