%% Batch Run StepWidth

% Influence of the kinematic constraint on stepwidth
clear all; close all; clc;

%% Default settings

% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 4;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% Folder with default functions
S.subject            = 's1_Poggensee';

% output folder
S.ResultsFolder     = 'Sens_HipF';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)
S.savename_ig   = 'NoExo';

% dataset with exoskeleton torque profile
S.DataSet       = 'PoggenSee2020_AFO';

%% Sensitivity weight metabolic rate
HipVect = 100:10:200;
nSim = length(HipVect);

for i =1:nSim   
    
    %scale factor for the force in the hip muscles
    Fhip = HipVect(i);
    
    % select the CasadiFolder
    S.CasadiFunc_Folders = ['Casadi_s1Pog_ScaleParam_HipF_' num2str(Fhip)];

    % normal walking simulation
    S.ExternalFunc  = 'PredSim_3D_Pog_s1_mtp.dll';        % this one is with the pinjoint mtp
    S.ExoBool       = 0;    S.ExoScale      = 0;    
    S.savename      = ['NoExo_HipF_' num2str(Fhip)];
    f_PredSim_PoggenSee2020(S);
    ct = ct+1;    

    % passive simulation
    S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % this one is with the pinjoint mtp
    S.ExoBool       = 1;    S.ExoScale      = 0;    
    S.savename      = ['Passive_HipF_' num2str(Fhip)];
    f_PredSim_PoggenSee2020(S);
    ct = ct+1;
    
    % active simulation
    S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % this one is with the pinjoint mtp
    S.ExoBool       = 1;    S.ExoScale      = 1;    
    S.savename      = ['Active_HipF_' num2str(Fhip)];
    f_PredSim_PoggenSee2020(S);
    ct = ct+1;
end