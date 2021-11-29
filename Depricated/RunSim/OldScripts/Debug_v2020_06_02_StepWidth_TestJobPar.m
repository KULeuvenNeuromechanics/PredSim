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
S.ResultsFolder     = 'Sens_StepWidth';

% select tendon stiffness of 20
S.CasadiFunc_Folders = 'Casadi_s1Pog_ScaleParam_k20';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)
S.savename_ig   = 'NoExo';

%% Open a cluster to run batch processes
% open cluster
myCluster = parcluster('Maarten_LocalProfile1');

% path information for cluster
StartPath = pwd;
MainPath = StartPath(1:end-7);
CasadiFiles = fullfile(MainPath,'CasADiFunctions',S.CasadiFunc_Folders);
PathPolynomials = fullfile(MainPath,'Polynomials',S.subject);
ExoPath = fullfile(MainPath,'Data','Poggensee_2020');
pathExternalFunctions = fullfile(MainPath,'ExternalFunctions');

%% Imposing change in stepwidth with passive exoskeleton

S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % this one is with the pinjoint mtp
S.ExoBool       = 1;
S.ExoScale      = 0;
S.DataSet       = 'PoggenSee2020_AFO';

SWidthV = 0:0.01:0.15;
nSim = length(SWidthV);

for i=1
    % set the kinematic constraint
    S.Constr.calcn = 0.09 + SWidthV(i);
    S.Constr.toes = 0.1 + SWidthV(i);
    S.Constr.tibia = 0.11 + SWidthV(i);
    % Change the output name
    WidthStr = round(SWidthV(i)*100);
    S.savename      = ['Passive_dFoot_' num2str(WidthStr) 'cm'];
    % run the simulation
    jobTest = batch(myCluster,'f_PredSim_PoggenSee2020_CalcnT',0,{S},...
        'CurrentFolder',StartPath,'AdditionalPaths',{CasadiFiles,PathPolynomials,ExoPath,pathExternalFunctions});
end

