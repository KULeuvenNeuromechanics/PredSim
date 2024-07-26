% Right and left muscles indices. Hard-coded for 3D model
% 'Falisse_et_al_2022.osim'. This may be automated
S.subject.idx_m_r = [1:43,87,89,91];
S.subject.idx_m_l = [44:86,88,90,92];

%% Example 1: Symmetric synergies, free weights
S.subject.synergies = 1; % 1 = implement muscle synergies
S.subject.TrackSynW = 0; % 1 = track synergy weights
S.misc.gaitmotion_type = 'HalfGaitCycle';
S.subject.NSyn = 4;

%% Example 2: Symmetric synergies, track weights 
% (from the case of 4 synergies in Pitto et al 2020)
S.subject.synergies = 1; % 1 = implement muscle synergies
S.subject.TrackSynW = 1; % 1 = track synergy weights
S.misc.gaitmotion_type = 'HalfGaitCycle';
S.subject.NSyn = 4; % NSyn_r = NSyn_l
S.subject.TrackSynW_NSyn_r = 4; % number of tracked synergies (may be different from the number of synergies)
S.subject.knownSynW_r = {'rect_fem_r', [0 0.28 0 0.75],...
    'vas_lat_r', [0.1 0.1 0.05 0.8],...
    {'bifemlh_r','bifemsh_r'}, [0.6 0.4 0.15 0.2],...
    {'semiten_r','semimem_r'}, [0 0 1 0],...
    'tib_ant_r', [0.2 1 0.01 0],...
    {'med_gas_r','lat_gas_r'}, [0.95 0.08 0.04 0.03],...
    'soleus_r', [1 0.05 0 0.05],...
    {'glut_med1_r','glut_med2_r','glut_med3_r'}, [0.15 0 0.04 1]};

%% Example 3: Non-symmetric synergies, free weights
S.subject.synergies = 1; % 1 = implement muscle synergies
S.subject.TrackSynW = 0; % 1 = track synergy weights
S.misc.gaitmotion_type = 'FullGaitCycle';
S.subject.NSyn_r = 5;
S.subject.NSyn_l = 2;

%% Example 4: Non-symmetric synergies, track weights for only left leg 
% (from experimental data from a CP child)
S.subject.synergies = 1; % 1 = implement muscle synergies
S.subject.TrackSynW = 1; % 1 = track synergy weights
S.misc.gaitmotion_type = 'FullGaitCycle';
S.subject.NSyn_r = 5;
S.subject.NSyn_l = 2;
S.subject.TrackSynW_side = 'onlyLeft'; % 'onlyRight'; % 'RightLeft';
S.subject.TrackSynW_NSyn_l = 2; % number of tracked synergies (may be different from the number of synergies)
S.subject.knownSynW_l = {'rect_fem_l', [0.3265 0.2298],...
    'vas_lat_l', [0.3214 0.1277],...
    {'bifemlh_l','bifemsh_l'}, [0.1940 0.3116],...
    {'semiten_l','semimem_l'}, [0.1393 0.3683],...
    'tib_ant_l', [0.0914 0.5599],...
    {'med_gas_l','lat_gas_l'}, [0.4883 0.4313],...
    'soleus_l', [0.6271 0.2091],...
    {'glut_med1_l','glut_med2_l','glut_med3_l'}, [0.2610 0.2617]};