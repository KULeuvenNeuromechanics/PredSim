
%% Batch Sim exoskeleton assistance

% Default simulations antoine with exoskeleton assistance.
% Note currently without inertia of the exoskeleton
clear all; close all; clc;

% add casadi to the matlab path (needed on simulation computer because
% default version is 3.4.5
rmpath(genpath('G:\PhD\Matlabcodes\Casadi\casadi'));
addpath(genpath('D:\MaartenAfschrift\softInstall\casadi'));
addpath(genpath('D:\MaartenAfschrift\GitProjects\3dpredictsim'));


%% Default settings

% flow control
S.Flow.solveProblem     = 0;   % set to 1 to solve problem
S.Flow.analyseResults   = 1;   % set to 1 to analyze results
S.Flow.loadResults      = 0;   % set to 1 to load results
S.Flow.saveResults      = 1;   % set to 1 to save sens. results
S.Flow.checkBoundsIG    = 0;   % set to 1 to visualize guess-bounds
S.Flow.writeIKmotion    = 1;   % set to 1 to write .mot file

% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.W.E       = 500;      % weight metabolic energy rate
S.W.Ak      = 50000;    % weight joint accelerations
S.W.ArmE    = 10^6;     % weight arm excitations
S.W.passMom = 1000;     % weight passive torques
S.W.A       = 2000;     % weight muscle activations
S.W.exp_E   = 2;        % power metabolic energy
S.W.Mtp     = 10^5;     % weight mtp excitations
S.W.u       = 0.001;    % weight on excitations arms actuators
S.IGsel     = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.ContactID = 1;        % contact model identifier
S.IGmodeID  = 1;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution)
S.IGcase    = 0;        % initial guess case identifier
S.h_weak    = 0;        % weakness hip actuators
S.Max_s     = 0;        % maximal contraction velocity identifier
S.pf_weak   = 0;        % weakness ankle plantaflexors
S.mE        = 0;        % metabolic energy model identifier
S.coCont    = 0;        % co-contraction identifier
S.NThreads  = 2;        % number of threads for parallel computing

% ipopt options
S.linear_solver = 'mumps';
S.tol_ipopt     = 4;

% quasi random initial guess, pelvis y position
% S.IG_PelvisY = 0.903;

% Point to the right subject
% S.CasadiFunc_Folders = 'Casadi_s1_Poggensee_k17mtp';
S.CasadiFunc_Folders = 'Casadi_subject1_k17mtp';
S.subject            = 'subject1';

% output folder
S.ResultsFolder = 'BatchSim_2020_03_10_DefaultSubj';

%% open cluster
if S.Flow.solveProblem
    myCluster = parcluster('Maarten_8Cores');
end


%% Passive exoskeleton
ct = 1; % counter 4 jobs
% external function
S.ExternalFunc  = 'PredSim_mtpPin_cm2.dll';        % this one is with the pinjoint mtp 
S.ExternalFunc2 = 'PredSim_mtpPin_pp_cm2.dll';    % this one is with the pinjoint mtp 
S.savename      = 'Passive';
S.loadname      = 'Passive';
S.ExoBool       = 1;
S.ExoScale      = 0;
S.DataSet       = 'PoggenSee2020_AFO';
S.IGmodeID      = 1;        % initial guess mode identifier
if S.Flow.solveProblem
    jobs(ct) = batch(myCluster,'f_PredSim_PoggenSee2020_DefaultS',0,{S}); ct =ct +1;
end

%% 100% assistance exoskeleton

% external function
S.ExternalFunc = 'PredSim_mtpPin_cm2.dll';        % this one is with the pinjoint mtp 
S.ExternalFunc2 = 'PredSim_mtpPin_pp_cm2.dll';    % this one is with the pinjoint mtp 

S.savename      = 'Active';
S.loadname      = 'Active';
S.ExoBool       = 1;
S.ExoScale      = 1;
S.DataSet       = 'PoggenSee2020_AFO';
S.IGmodeID      = 1;        % initial guess mode identifier
if S.Flow.solveProblem
    jobs(ct) = batch(myCluster,'f_PredSim_PoggenSee2020_DefaultS',0,{S});ct =ct +1;
end
%% Loop over level of assistance
% external function
S.ExternalFunc  = 'PredSim_mtpPin_cm2.dll';        % this one is with the pinjoint mtp 
S.ExternalFunc2 = 'PredSim_mtpPin_pp_cm2.dll';    % this one is with the pinjoint mtp 
AVect = 0:20:200;
for j = 1:length(AVect)
    S.savename      = ['Active_' num2str(AVect(j))];
    S.loadname      = ['Active_' num2str(AVect(j))];
    S.ExoBool       = 1;
    S.ExoScale      = AVect(j)./100;
    S.DataSet       = 'PoggenSee2020_AFO';
    S.IGmodeID      = 1;        % initial guess mode identifier
    if S.Flow.solveProblem
        jobs(ct) = batch(myCluster,'f_PredSim_PoggenSee2020_DefaultS',0,{S});ct =ct +1;
    end
end

%% Vector with names
RFolder = S.ResultsFolder;
dPath = fullfile('C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\Results',RFolder);

% Create Vector with Names
Names = {'Passive','Active'};
for j = 1:length(AVect)
   Names{j+2} =  ['Active_' num2str(AVect(j))];
end

%% Post processing

for i = 1:length(Names)
    if exist(fullfile(dPath,[Names{i} '.mat']),'file')
        f_LoadSim_PoggenSee2020_DefaultS(RFolder,Names{i});
    end
end

%% Plot Results

% % No Exo in black
% Cs = [0 0 0];
% PlotResults_3DSim(fullfile(dPath,[Names{1} '_pp.mat']),Cs);
% h = gcf;

% passive in red
Cs = [1 0 0];
PlotResults_3DSim(fullfile(dPath,[Names{1} '_pp.mat']),Cs,h);

% active in blue
Cs = [0 0 1];
PlotResults_3DSim(fullfile(dPath,[Names{2} '_pp.mat']),Cs,h);

% percentage of assistanc ein copper
CsVect = copper(length(AVect));
for i=1:length(AVect)
    Cs = CsVect(i,:);
    PlotResults_3DSim(fullfile(dPath,[Names{i+2} '_pp.mat']),Cs,h);
end
