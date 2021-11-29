%% Batch Run StepWidth

% Influence of the kinematic constraint on stepwidth
clear all; close all; clc;

%% Default settings

% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 8;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% Folder with default functions
S.subject            = 's1_Poggensee';

% output folder
S.ResultsFolder     = 'Simulation_Browning2008';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)
S.savename_ig   = 'NoExo';

% dataset with exoskeleton torque profile
S.DataSet       = 'PoggenSee2020_AFO';

% select the CasadiFolder
S.CasadiFunc_Folders = 'Casadi_s1Pog_mtp';

% add post processing external function (makes it a lot easier for
% postprocessing)
S.ExternalFunc2 = 'Browning_2008_pp.ddl';   

%% Run with all external functions
ExtFuncs = {'Browning_2008_Reference','Browning_2008_Femur2','Browning_2008_Femur4','Browning_2008_Femur8',...
    'Browning_2008_Foot2','Browning_2008_Foot4','Browning_2008_Pelvis4','Browning_2008_Pelvis8','Browning_2008_Pelvis12',...
    'Browning_2008_Pelvis16','Browning_2008_Tibia2','Browning_2008_Tibia4'};
nSim = length(ExtFuncs);
for i =1:nSim
    % normal walking simulation
    FuncSel = ExtFuncs{i};
    S.ExternalFunc  = [FuncSel '.dll'];        % this one is with the pinjoint mtp
    S.ExoBool       = 0;    
    S.ExoScale      = 0;
    S.savename      = FuncSel;
    f_PredSim_PoggenSee2020(S);
end

