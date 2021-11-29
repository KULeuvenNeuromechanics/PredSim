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
S.ResultsFolder     = 'Sensitivity_ContactSpheres';

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

% set model name
S.ModelName = 'Gait92';

%% Run with all external functions

% external functions
ExtFuncs = {'SimExo_3D_talus_out','SimExo_3D_c1','SimExo_3D_c2','SimExo_3D_c5'};
ExtFunc2 = {'SimExo_3D_ExportAll','SimExo_3D_c1_pp','SimExo_3D_c2_pp','SimExo_3D_c5_pp'};

for i=1:length(ExtFuncs)
    % normal walking simulation
    FuncSel         = ExtFuncs{i};    
    S.ExternalFunc  = [FuncSel '.dll'];        % this one is with the pinjoint mtp
    S.ExternalFunc2 = [ExtFunc2{i} '.dll'];
    
    % passive exoskeleton
    S.ExoBool       = 0;
    S.ExoScale      = 0;
    S.savename      = ['P_ ' FuncSel];
    f_PredSim_Gait92(S);

    % active exoskeleton
    S.ExoBool       = 1;
    S.ExoScale      = 1;
    S.savename      = ['A_ ' FuncSel];
    f_PredSim_Gait92(S);
end



