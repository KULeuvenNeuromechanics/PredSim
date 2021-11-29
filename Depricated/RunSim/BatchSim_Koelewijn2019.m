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

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% Folder with default functions
S.subject            = 's1_Poggensee';

% output folder
S.ResultsFolder     = 'Simulation_Koelewijn2019';

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
S.mass = 64;    % mass in kg

% add post processing external function (makes it a lot easier for
% postprocessing)
S.ExternalFunc2 = 'Analyse_s1Pog.ddl';

% set model name
S.ModelName = 'Gait92';

%% Run with all external functions
ExtFuncs = {'Flat_s1Pog','Decl4_s1Pog','Decl8_s1Pog','Incl4_s1Pog','Incl8_s1Pog'};
Speeds = [0.8 1.3];
nSim = length(ExtFuncs);
for s = 1:2
    speedSel = Speeds(s);
    speedName = round(speedSel*10);
    for i =1:nSim
        % normal walking simulation
        FuncSel = ExtFuncs{i};
        S.ExternalFunc  = [FuncSel '.dll'];        % this one is with the pinjoint mtp
        S.ExoBool       = 0;
        S.ExoScale      = 0;
        S.savename      = [FuncSel '_speed_' num2str(speedName)];
        S.v_tgt     	= speedSel;
        f_PredSim_Gait92(S);
    end
end



