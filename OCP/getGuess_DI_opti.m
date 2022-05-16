function guess = getGuess_DI_opti(S,model_info,scaling,d)
% --------------------------------------------------------------------------
% getGuess_DI_opti
%   This script provides an inital guess for the design variables.
%   The guess is data-informed (DI). We use experimental data to provide an
%   initial guess of the joint variables but set constant values to the 
%   muscle variable and the arm variables. We use a pre-defined final time 
%   that is function of the imposed speed.
%   
% INPUT:
%   - S -
%   * setting structure S
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
%   - scaling -
%   * 
% 
%   - d -
%   * 
%
% OUTPUT:
%   - guess -
%   * 
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

N = S.solver.N_meshes; % number of mesh intervals
nq = model_info.ExtFunIO.jointi.nq;
NMuscle = model_info.muscle_info.NMuscle;
coordinate_names = model_info.ExtFunIO.coord_names.all;

%%

Qs_guess = getIK(S.subject.IG_selection,model_info);
if strcmp(S.misc.gaitmotion_type,'FullGaitCycle')
    endIdx = round(size(Qs_guess.allfilt,1)*100/S.subject.IG_selection_gaitCyclePercent);
elseif strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    endIdx = round(size(Qs_guess.allfilt,1)*50/S.subject.IG_selection_gaitCyclePercent);
end
Qs_guess_IG.allfilt = Qs_guess.allfilt(1:endIdx,:);
Qs_guess_IG.time = Qs_guess.time(1:endIdx,:);
Qs_guess_IG.colheaders = Qs_guess.colheaders;
time_IC = [Qs_guess_IG.time(1),Qs_guess_IG.time(end)];

%% Spline approximation of Qs to get Qdots and Qdotdots
Qs_spline.data = zeros(size(Qs_guess_IG.allfilt));
Qs_spline.data(:,1) = Qs_guess_IG.allfilt(:,1);
Qdots_spline.data = zeros(size(Qs_guess_IG.allfilt));
Qdots_spline.data(:,1) = Qs_guess_IG.allfilt(:,1);
Qdotdots_spline.data = zeros(size(Qs_guess_IG.allfilt));
Qdotdots_spline.data(:,1) = Qs_guess_IG.allfilt(:,1);
for i = 2:size(Qs_guess_IG.allfilt,2)
    Qs_guess_IG.datafiltspline(i) = spline(Qs_guess_IG.allfilt(:,1),Qs_guess_IG.allfilt(:,i));
    [Qs_spline.data(:,i),Qdots_spline.data(:,i),...
        Qdotdots_spline.data(:,i)] = ...
        SplineEval_ppuval(Qs_guess_IG.datafiltspline(i),Qs_guess_IG.allfilt(:,1),1);
end

for i=1:nq.all
    coordinate = coordinate_names{i};
    coord_idx = model_info.ExtFunIO.coordi.(coordinate);
    guess.Qs_all(:,coord_idx) = Qs_spline.data(:,strcmp(Qs_guess_IG.colheaders(1,:),coordinate));
    guess.Qdots_all(:,coord_idx) = Qdots_spline.data(:,strcmp(Qs_guess_IG.colheaders(1,:),coordinate));
    guess.Qdotdots_all(:,coord_idx) = Qdotdots_spline.data(:,strcmp(Qs_guess_IG.colheaders(1,:),coordinate));
end
    
% Interpolation
Qs_time = Qs_spline.data(:,strcmp(Qs_guess_IG.colheaders(1,:),'time'));
time_expi.Qs(1) = find(round(Qs_time,3) == round(time_IC(1),3));
time_expi.Qs(2) = find(round(Qs_time,3) == round(time_IC(2),3));
step = (Qs_time(time_expi.Qs(2))-Qs_time(time_expi.Qs(1)))/(N-1);
interval = Qs_time(time_expi.Qs(1)):step:Qs_time(time_expi.Qs(2));
guess.Qs = interp1(round(Qs_time,4),guess.Qs_all,round(interval,4));
guess.Qs(:,model_info.ExtFunIO.jointi.base_forward) = guess.Qs(:,model_info.ExtFunIO.jointi.base_forward) - ....
    guess.Qs(1,model_info.ExtFunIO.jointi.base_forward);

if S.subject.adapt_IG_pelvis_y
    % Adjust pelvis height
    guess.Qs(:,model_info.ExtFunIO.coordi.pelvis_ty) = guess.Qs(:,model_info.ExtFunIO.coordi.pelvis_ty) ...
        - mean(guess.Qs(:,model_info.ExtFunIO.coordi.pelvis_ty)) + model_info.IG_pelvis_y;
end

% Interpolation
guess.Qdots = interp1(round(Qs_time,4),guess.Qdots_all,round(interval,4));

% Interpolation
guess.Qdotdots = interp1(round(Qs_time,4),guess.Qdotdots_all,round(interval,4));

%% Muscle variables
guess.a = 0.1*ones(N,NMuscle);
guess.vA = 0.01*ones(N,NMuscle);
guess.FTtilde = 0.1*ones(N,NMuscle);
guess.dFTtilde = 0.01*ones(N,NMuscle);

%% Arm activations
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
    dx = guess.Qs(end,model_info.ExtFunIO.jointi.base_forward) - ...
        guess.Qs(end-1,model_info.ExtFunIO.jointi.base_forward);
    inv_X_Qs(model_info.ExtFunIO.jointi.base_forward) = ...
        guess.Qs(end,model_info.ExtFunIO.jointi.base_forward) + dx;

    guess.Qs = [guess.Qs; inv_X_Qs];
    guess.Qdots = [guess.Qdots; inv_X_Qdots];
    guess.a = [guess.a; guess.a(1,model_info.ExtFunIO.symQs.MusInvB)];
    guess.FTtilde = [guess.FTtilde; guess.FTtilde(1,model_info.ExtFunIO.symQs.MusInvB)];
    guess.a_a = [guess.a_a; guess.a_a(1,:)];
else
    dx = guess.Qs(end,model_info.ExtFunIO.jointi.base_forward) - ...
        guess.Qs(end-1,model_info.ExtFunIO.jointi.base_forward);
    guess.Qs = [guess.Qs; guess.Qs(1,:)];
    guess.Qs(end,model_info.ExtFunIO.jointi.base_forward) = guess.Qs(end-1,model_info.ExtFunIO.jointi.base_forward) + dx;
    guess.Qdots = [guess.Qdots; guess.Qdots(1,:)];
    guess.a = [guess.a; guess.a(1,:)];
    guess.FTtilde = [guess.FTtilde; guess.FTtilde(1,:)];
    guess.a_a = [guess.a_a; guess.a_a(1,:)];
end


%% Final time
% The final time is function of the imposed speed
all_speeds = 0.73:0.1:5;
% all_speeds = 0.73:0.1:2.73;
all_tf = 0.70:-((0.70-0.35)/(length(all_speeds)-1)):0.35;
idx_speed = find(all_speeds==S.subject.v_pelvis_x_trgt);
if isempty(idx_speed)
    idx_speed = find(all_speeds > S.subject.v_pelvis_x_trgt,1,'first');
end
guess.tf = all_tf(idx_speed);

%% Scaling
guess.Qs = guess.Qs./repmat(scaling.Qs,N+1,1);
guess.Qdots = guess.Qdots./repmat(scaling.Qdots,N+1,1);
guess.Qdotdots = guess.Qdotdots./repmat(scaling.Qdotdots,N,1);
guess.a         = (guess.a)./repmat(scaling.a,N+1,size(guess.a,2));
guess.FTtilde   = (guess.FTtilde)./repmat(scaling.FTtilde,N+1,1);
guess.vA        = (guess.vA)./repmat(scaling.vA,N,size(guess.vA,2));
guess.dFTtilde  = (guess.dFTtilde)./repmat(scaling.dFTtilde,N,...
    size(guess.dFTtilde,2));


%% Collocation points
guess.a_col = zeros(d*N,NMuscle);
guess.FTtilde_col = zeros(d*N,NMuscle);
guess.Qs_col = zeros(d*N,nq.all);
guess.Qdots_col = zeros(d*N,nq.all);
guess.a_a_col = zeros(d*N,nq.torqAct);
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
end
