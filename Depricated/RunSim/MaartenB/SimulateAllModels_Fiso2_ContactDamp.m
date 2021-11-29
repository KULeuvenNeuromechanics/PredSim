clear all; close all; clc;

% settings for optimization
S.v_tgt     = 4;     % average speed
S.N         = 30;       % number of mesh intervals
S.NThreads  = 8;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% output folder
S.ResultsFolder     = 'AllModels_MaartenB_ContactDamp_V2';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 2;        % initial guess mode identifier (1 walk, 2 run, 3 prev.solution, 4 solution from /IG/Data folder)
S.IKfile_guess = 'OpenSimModel\IK_Guess_Default.mat'; % or OpenSimModel\IK_Guess_Default_Running.mat
S.IKfile_Bounds = 'OpenSimModel\IK_Bounds_Default.mat'; % or OpenSimModel\IK_Guess_Default_Running.mat

% Select model
S.ModelName = 'Gait92'; % other option is 'Rajagopal'

% mass of the subject
S.mass = 64;    % mass in kg (needed ?)

% bounds based on running
S.Bounds_Running = true;

% loop over the 6 models
OutNames = {'Sedgeford_Large','Sedgeford_Mean','Sedgeford_Small','Haya_Large','Haya_Mean','Haya_Small'};
OutNames_Fiso = {'Sedgeford_Large_Fiso2','Sedgeford_Mean_Fiso2',...
    'Sedgeford_Small_Fiso2','Haya_Large_Fiso2','Haya_Mean_Fiso2','Haya_Small_Fiso2'};
ExtFuncs = {'s1Pog_dis2','s1Pog_dis1','s1Pog_dis05','s1Pog_dis02'};

for j = 1:length(ExtFuncs)
    for i=1%:length(OutNames)
        FuncSel         = ExtFuncs{j};
        S.ExternalFunc  = [FuncSel '.dll'];
        S.ExternalFunc2  = [FuncSel '_pp.dll'];
        % selection folder with Casadi Functions
        S.CasadiFunc_Folders = OutNames_Fiso{i};
        % select folder with polynomials
        S.PolyFolder = OutNames{i};
        % savename
        S.savename   = FuncSel;
        % run the optimization
        f_PredSim_Gait92(S);     % run the optimization
%         f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-proces simulation results
    end
end
