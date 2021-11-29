%% Batch Run StepWidth

% reproducing experimental work of Anne Koelewijn in 2019:
% Metabolic cost calculations of gait using musculoskeletal energy models, a
% comparison study

% Influence of the kinematic constraint on stepwidth
clear all; close all; clc;

%% Default settings

% settings for optimization
% S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 8;        % number of threads for parallel computing

% walking speed
S.v_tgt     	= 1.25;

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% Folder with default functions
S.subject            = 's1_Poggensee';

% output folder
S.ResultsFolder     = 'Sensitivity_DisPCoeff';

% select the CasadiFolder
S.CasadiFunc_Folders = 'Casadi_s1Pog_mtp';

% select folder with polynomials
S.PolyFolder = 's1_Poggensee';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)
S.savename_ig   = 'NoExo';

% dataset with exoskeleton torque profile
S.DataSet       = 'PoggenSee2020_AFO';

% Select model
S.ModelName = 'Gait92'; % other option is 'Rajagopal'

% mass of the subject
S.mass = 68;    % mass in kg

% add post processing external function (makes it a lot easier for
% postprocessing)
S.ExternalFunc2 = 'Analyse_s1Pog.ddl';

% set model name
S.ModelName = 'Gait92';

%% Run with all external functions

% external functions
ExtFuncs = {'s1Pog_dis2','s1Pog_dis1','s1Pog_dis3','s1Pog_dis05','s1Pog_dis02','s1Pog_dis01','s1Pog_dis4'};
for i=1:length(ExtFuncs)
    % normal walking simulation
    FuncSel         = ExtFuncs{i};
    S.ExternalFunc  = [FuncSel '.dll'];        % this one is with the pinjoint mtp
    S.ExoBool       = 0;
    S.ExoScale      = 0;
    S.savename      = FuncSel;
    f_PredSim_Gait92(S);
end



