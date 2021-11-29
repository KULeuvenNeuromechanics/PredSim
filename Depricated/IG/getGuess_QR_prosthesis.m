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
function guess = getGuess_QR_prosthesis(N,nq,NMuscle,scaling,v_tgt,jointi)

%% Final time
% The final time is function of the imposed speed
all_speeds = 0.73:0.1:2.73;
all_tf = 0.70:-((0.70-0.35)/(length(all_speeds)-1)):0.35;
all_tf = all_tf*2;
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
guess.Qs(:,jointi.pelvis.ty) = 0.9385;

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

%% Scaling
guess.QsQdots   = guess.QsQdots./repmat(scaling.QsQdots,N,1);
guess.Qdotdots  = guess.Qdotdots./repmat(scaling.Qdotdots,N,1);
guess.a         = (guess.a)./repmat(scaling.a,N,size(guess.a,2));
guess.FTtilde   = (guess.FTtilde)./repmat(scaling.FTtilde,N,1);
guess.vA        = (guess.vA)./repmat(scaling.vA,N,size(guess.vA,2));
guess.dFTtilde  = (guess.dFTtilde)./repmat(scaling.dFTtilde,N,...
    size(guess.dFTtilde,2));

end
