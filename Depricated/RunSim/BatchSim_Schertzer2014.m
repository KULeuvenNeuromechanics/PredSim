%% Batch Run StepWidth

% Influence of the kinematic constraint on stepwidth
clear all; close all; clc;

%% Default settings

% settings for optimization
S.N         = 50;       % number of mesh intervals
S.NThreads  = 4;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% Folder with default functions
S.subject            = 's1_Poggensee';

% select folder with polynomials
S.PolyFolder = 's1_Poggensee';

% output folder
S.ResultsFolder     = 'Simulation_Schertzer';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)
S.savename_ig   = 'NoExo';

% dataset with exoskeleton torque profile
S.DataSet       = 'PoggenSee2020_AFO';

% mass of the subject
S.mass = 64;    % mass in kg

% select the CasadiFolder
S.CasadiFunc_Folders = 'Casadi_s1Pog_mtp';

% add post processing external function (makes it a lot easier for
% postprocessing)
S.ExternalFunc2 = 'Schertzer2014_pp.ddl';

% Select model
S.ModelName = 'Gait92'; % other option is 'Rajagopal'


%% Run with all external functions
ExtFuncs = {'Schertzer2014_Reference','Schertzer2014_Ankle1','Schertzer2014_Ankle2','Schertzer2014_Ankle05',...
    'Schertzer2014_Back2','Schertzer2014_Back7','Schertzer2014_Back10','Schertzer2014_Back16','Schertzer2014_Back22',...
    'Schertzer2014_Knee1','Schertzer2014_Knee2','Schertzer2014_Knee05'};
nSim = length(ExtFuncs);
Speedskmh = [4 5 6];
nSpeeds = length(Speedskmh);
% for j=1:nSpeeds
%     for i =1:nSim
%         % simulation settings
%         S.v_tgt     = Speedskmh(j)/3.6;     % average speed
%         FuncSel = ExtFuncs{i};
%         S.ExternalFunc  = [FuncSel '.dll'];        % this one is with the pinjoint mtp
%         S.ExoBool       = 0;
%         S.ExoScale      = 0;
%         S.savename      = [FuncSel '_' num2str(Speedskmh(j)) 'km_h'];
%         f_PredSim_Gait92(S);
%     end
% end
j = 1;
i = 1;
S.v_tgt     = Speedskmh(j)/3.6;     % average speed
FuncSel = ExtFuncs{i};
S.ExternalFunc  = [FuncSel '.dll'];        % this one is with the pinjoint mtp
S.ExoBool       = 0;
S.ExoScale      = 0;
S.savename      = [FuncSel '_' num2str(Speedskmh(j)) 'km_h'];
f_PredSim_Gait92(S);
