% Right and left muscles indices. Hard-coded for 3D model
% 'Falisse_et_al_2022'. This may be automated
S.subject.idx_m_r = [1:43,87,89,91];
S.subject.idx_m_l = [44:86,88,90,92];

%% Example 1: Symmetric synergies, free weights
S.subject.synergies = 1; % 1 = implement muscle synergies
S.subject.TrackSynW = 0; % 1 = track synergy weights
S.misc.gaitmotion_type = 'HalfGaitCycle';
S.subject.NSyn = 4; 

%% Example 2: Symmetric synergies, track weights (from Pitto et al 2020)
% The specific weights to be tracked and the muscle indices from the known
% weights is for now hard-coded here
S.subject.synergies = 1; % 1 = implement muscle synergies
S.subject.TrackSynW = 1; % 1 = track synergy weights
S.misc.gaitmotion_type = 'HalfGaitCycle';
S.subject.NSyn = 4; % NSyn_r = NSyn_l
Pitto2020_4Syn = [0 0.1 0.6 0 0.2 0.95 1 0.15;...
    0.28 0.1 0.4 0 1 0.08 0.05 0;...
    0 0.05 0.15 1 0.01 0.04 0 0.04;...
    0.75 0.8 0.2 0 0 0.03 0.05 1];
S.subject.knownSynW = [Pitto2020_4Syn(:,1),...
    Pitto2020_4Syn(:,2),...
    Pitto2020_4Syn(:,3),Pitto2020_4Syn(:,3),...
    Pitto2020_4Syn(:,4),Pitto2020_4Syn(:,4),...
    Pitto2020_4Syn(:,5),...
    Pitto2020_4Syn(:,6),Pitto2020_4Syn(:,6),...
    Pitto2020_4Syn(:,7),...
    Pitto2020_4Syn(:,8),Pitto2020_4Syn(:,8),Pitto2020_4Syn(:,8)]';
S.subject.knownSynW_idx = [28 31 9 10 7 8 38 32 33 34 1 2 3];
% 1 Rectus femoris,
% 2 Vastus lateralis,
% 3 Biceps femoris, (long head and short head)
% 4 medial hamstrings,  (semiten and semimem)
% 5 Tibialis anterior,
% 6 Gastrocnemius, (lat and med)
% 7 Soleus,
% 8 Gluteus medius (1,2,3)

%% Example 3: Non-symmetric synergies, free weights
S.subject.synergies = 1; % 1 = implement muscle synergies
S.subject.TrackSynW = 0; % 1 = track synergy weights
S.misc.gaitmotion_type = 'FullGaitCycle';
S.subject.NSyn_r = 5;
S.subject.NSyn_l = 2;

%% Example 4: Non-symmetric synergies, track weights for only left leg (from
% experimental data from a CP child)
% The specific weights to be tracked and the muscle indices from the known
% weights is for now hard-coded here
S.subject.synergies = 1; % 1 = implement muscle synergies
S.subject.TrackSynW = 1; % 1 = track synergy weights
S.misc.gaitmotion_type = 'FullGaitCycle';
S.subject.NSyn_r = 5;
S.subject.NSyn_l = 2;
S.subject.TrackSynW_side = 'onlyLeft'; % 'onlyRight'; % 'RightLeft';
CP_SynW = [0.3265    0.3214    0.1940    0.1393    0.0914    0.4883    0.6271    0.2610; ...
    0.2298    0.1277    0.3116    0.3683    0.5599    0.4313    0.2091    0.2617]; % from a specific CP child
S.subject.knownSynW_l = [CP_SynW(:,1),...
    CP_SynW(:,2),...
    CP_SynW(:,3),CP_SynW(:,3),...
    CP_SynW(:,4),CP_SynW(:,4),...
    CP_SynW(:,5),...
    CP_SynW(:,6),CP_SynW(:,6),...
    CP_SynW(:,7),...
    CP_SynW(:,8),CP_SynW(:,8),CP_SynW(:,8)]';
S.subject.knownSynW_idx_l = [28 31 9 10 7 8 38 32 33 34 1 2 3]; % indeces start at 1 for the left muscles,
% as well as for the right muscles
% 1 Rectus femoris,
% 2 Vastus lateralis,
% 3 Biceps femoris, (long head and short head)
% 4 medial hamstrings,  (semiten and semimem)
% 5 Tibialis anterior,
% 6 Gastrocnemius, (lat and med)
% 7 Soleus,
% 8 Gluteus medius (1,2,3)