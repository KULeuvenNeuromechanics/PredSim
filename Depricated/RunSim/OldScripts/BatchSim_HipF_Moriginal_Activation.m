%% Batch Run StepWidth

% Influence of the kinematic constraint on stepwidth
clear all; close all; clc;

%% Default settings

% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 2;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% Folder with default functions
S.subject            = 's1_Poggensee';

% output folder
S.ResultsFolder     = 'Sens_HipF_mOr_act';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)
S.savename_ig   = 'NoExo';

% dataset with exoskeleton torque profile
S.DataSet       = 'PoggenSee2020_AFO';

%% Open a cluster to run batch processes
% open cluster
myCluster = parcluster('Maarten_LocalProfile1');


%% Sensitivity weight metabolic rate
% scale factor for hip muscles
HipScaleFactor = 0.6:0.2:3; % increase strength of hip muscles
nSim = length(HipScaleFactor);

ct = 1;
for i =1:nSim   
    
    
    %scale factor for the force in the hip muscles
    HipScaleSel = HipScaleFactor(i);
    HipScaleSelSTR = round(HipScaleSel*100);
    
    % select the CasadiFolder
    S.CasadiFunc_Folders = ['Casadi_s1Pog_HipF_Morig_' num2str(HipScaleSelSTR)];
    
    % baseline activation muscles
    S.Bounds.ActLower = 0.05;           % set baseline activity to 0.02
    S.Bounds.ActLowerHip = 0.05./HipScaleSel;    % only for the hip joint
    
    % information for cluster
    StartPath = pwd;
    MainPath = StartPath(1:end-7);
    CasadiFiles = fullfile(MainPath,'CasADiFunctions',S.CasadiFunc_Folders);
    PathPolynomials = fullfile(MainPath,'Polynomials',S.subject);
    ExoPath = fullfile(MainPath,'Data','Poggensee_2020');
    pathExternalFunctions = fullfile(MainPath,'ExternalFunctions');    

    % normal walking simulation
    S.ExternalFunc  = 'PredSim_3D_Pog_s1_mtp.dll';        % this one is with the pinjoint mtp
    S.ExoBool       = 0;    S.ExoScale      = 0;    
    S.savename      = ['NoExo_HipF_' num2str(HipScaleSelSTR)];
    jobs(ct) = batch(myCluster,'f_PredSim_PoggenSee2020',0,{S},...
       'CurrentFolder',StartPath,'AdditionalPaths',{CasadiFiles,PathPolynomials,ExoPath,pathExternalFunctions});
    ct = ct+1;%     f_PredSim_PoggenSee2020(S);
    
    % passive simulation
    S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % this one is with the pinjoint mtp
    S.ExoBool       = 1;    S.ExoScale      = 0;    
    S.savename      = ['Passive_HipF_' num2str(HipScaleSelSTR)];

    jobs(ct) = batch(myCluster,'f_PredSim_PoggenSee2020',0,{S},...
       'CurrentFolder',StartPath,'AdditionalPaths',{CasadiFiles,PathPolynomials,ExoPath,pathExternalFunctions});
    ct = ct+1;%     f_PredSim_PoggenSee2020(S);
    
    % active simulation
    S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % this one is with the pinjoint mtp
    S.ExoBool       = 1;    S.ExoScale      = 1;    
    S.savename      = ['Active_HipF_' num2str(HipScaleSelSTR)];
    jobs(ct) = batch(myCluster,'f_PredSim_PoggenSee2020',0,{S},...
       'CurrentFolder',StartPath,'AdditionalPaths',{CasadiFiles,PathPolynomials,ExoPath,pathExternalFunctions});
    ct = ct+1;    %     f_PredSim_PoggenSee2020(S);

end