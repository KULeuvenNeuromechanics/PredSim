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
S.ResultsFolder     = 'Sens_Obj';

% select tendon stiffness of 20
S.CasadiFunc_Folders = 'Casadi_s1Pog_ScaleParam_k20';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)
S.savename_ig   = 'NoExo';

% external function
S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % this one is with the pinjoint mtp

% dataset with exoskeleton torque profile
S.DataSet       = 'PoggenSee2020_AFO';

%% Open a cluster to run batch processes
% open cluster
% myCluster = parcluster('Maarten_LocalProfile1');
% 
% % path information for cluster
% StartPath = pwd;
% MainPath = StartPath(1:end-7);
% CasadiFiles = fullfile(MainPath,'CasADiFunctions',S.CasadiFunc_Folders);
% PathPolynomials = fullfile(MainPath,'Polynomials',S.subject);
% ExoPath = fullfile(MainPath,'Data','Poggensee_2020');
% pathExternalFunctions = fullfile(MainPath,'ExternalFunctions');

%% Sensitivity weight metabolic rate
% ct = 1;
% EVect = [0 10 20 50 100 200 300 500 1000 2000 5000];
% nSim = length(EVect);
% S.W.Ak      = 50000;    % weight joint accelerations
% S.W.A       = 2000;     % weight muscle activations
% for i =1:nSim   
%     % set weight on joint accelerations
%     S.W.E      = EVect(i);    % weight joint accelerations
%     % passive simulations
%     S.ExoBool       = 1;    S.ExoScale      = 0;    
%     S.savename      = ['Passive_metab_' num2str(S.W.Ak)];
%     f_PredSim_PoggenSee2020_CalcnT(S);
%     %jobs(ct) = batch(myCluster,'f_PredSim_PoggenSee2020_CalcnT',0,{S},...
%     %    'CurrentFolder',StartPath,'AdditionalPaths',{CasadiFiles,PathPolynomials,ExoPath,pathExternalFunctions});
%     ct = ct+1;
%     % active simulation
%     S.ExoBool       = 1;    S.ExoScale      = 1;    
%     S.savename      = ['Active_metab_' num2str(S.W.Ak)];
%     f_PredSim_PoggenSee2020_CalcnT(S);
% %     jobs(ct) = batch(myCluster,'f_PredSim_PoggenSee2020_CalcnT',0,{S},...
% %         'CurrentFolder',StartPath,'AdditionalPaths',{CasadiFiles,PathPolynomials,ExoPath,pathExternalFunctions});
%     ct = ct+1;
% end
%% Sensitivity weight joint accelerations
ct = 1;
QddVect = [0 50 100 1000 10000 20000 50000 80000 100000 200000];
nSim = length(QddVect);
S.W.E       = 500;      % weight metabolic energy rate
S.W.A       = 2000;     % weight muscle activations
for i =1:nSim   
    % set weight on joint accelerations
    S.W.Ak      = QddVect(i);    % weight joint accelerations
    % passive simulations
    S.ExoBool       = 1;    S.ExoScale      = 0;    
    S.savename      = ['Passive_qdd_' num2str(S.W.Ak)];
    f_PredSim_PoggenSee2020_CalcnT(S);
    %jobs(ct) = batch(myCluster,'f_PredSim_PoggenSee2020_CalcnT',0,{S},...
    %    'CurrentFolder',StartPath,'AdditionalPaths',{CasadiFiles,PathPolynomials,ExoPath,pathExternalFunctions});
    ct = ct+1;
    % active simulation
    S.ExoBool       = 1;    S.ExoScale      = 1;    
    S.savename      = ['Active_qdd_' num2str(S.W.Ak)];
    f_PredSim_PoggenSee2020_CalcnT(S);
    %jobs(ct) = batch(myCluster,'f_PredSim_PoggenSee2020_CalcnT',0,{S},...
    %    'CurrentFolder',StartPath,'AdditionalPaths',{CasadiFiles,PathPolynomials,ExoPath,pathExternalFunctions});
    ct = ct+1;
end

%% Sensitivity weight activations

% ct = 1;
% AVect = [0 10 100 200 500 1000 2000 5000 10000 20000];
% nSim = length(AVect);
% S.W.E       = 500;      % weight metabolic energy rate
% S.W.Ak      = 50000;    % weight joint accelerations
% for i =1:nSim   
%     % set weight on joint accelerations
%     S.W.A      = AVect(i);    % weight joint accelerations
%     % passive simulations
%     S.ExoBool       = 1;    S.ExoScale      = 0;    
%     S.savename      = ['Passive_act_' num2str(S.W.A)];
%     f_PredSim_PoggenSee2020_CalcnT(S);
%     %jobs(ct) = batch(myCluster,'f_PredSim_PoggenSee2020_CalcnT',0,{S},...
%     %    'CurrentFolder',StartPath,'AdditionalPaths',{CasadiFiles,PathPolynomials,ExoPath,pathExternalFunctions});
%     ct = ct+1;
%     % active simulation
%     S.ExoBool       = 1;    S.ExoScale      = 1;    
%     S.savename      = ['Active_act_' num2str(S.W.A)];
%     f_PredSim_PoggenSee2020_CalcnT(S);
%     %jobs(ct) = batch(myCluster,'f_PredSim_PoggenSee2020_CalcnT',0,{S},...
%     %    'CurrentFolder',StartPath,'AdditionalPaths',{CasadiFiles,PathPolynomials,ExoPath,pathExternalFunctions});
%     ct = ct+1;
% end



