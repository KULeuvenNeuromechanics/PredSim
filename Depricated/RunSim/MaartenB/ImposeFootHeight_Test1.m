clear all; close all; clc;

% settings for optimization
S.v_tgt     = 4;     % average speed
S.N         = 30;       % number of mesh intervals
S.NThreads  = 4;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% output folder
S.ResultsFolder     = 'AllModels_MaartenB_FootHeight';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
% S.IGmodeID      = 2;        % initial guess mode identifier (1 walk, 2 run, 3 prev.solution, 4 solution from /IG/Data folder)
% S.IKfile_guess = 'OpenSimModel\IK_Guess_Default.mat'; % or OpenSimModel\IK_Guess_Default_Running.mat
% S.IKfile_Bounds = 'OpenSimModel\IK_Bounds_Default.mat'; % or OpenSimModel\IK_Guess_Default_Running.mat

S.IGmodeID      = 3;        % initial guess mode identifier (1 walk, 2 run, 3 prev.solution, 4 solution from /IG/Data folder)
S.ResultsF_ig   = 'AllModels_MaartenB_FootHeight';
S.savename_ig   = 'Sedgeford_Mean_tf_FootHeight';

% Select model
S.ModelName = 'Gait92'; % other option is 'Rajagopal'

% mass of the subject
S.mass = 64;    % mass in kg (needed ?)

% Simulation without exoskeelton
S.ExternalFunc  = 'S1Pog_Running_V2.dll';        % external function
S.ExternalFunc2 = 'S1Pog_Running_pp.dll';     % external function for post-processing

% bounds based on running
S.Bounds_Running = true;

% impose stride frequency
FinalTime = 0.35;

% selection folder with Casadi Functions
S.CasadiFunc_Folders = 'Sedgeford_Mean_Fiso2';

% select folder with polynomials
S.PolyFolder = 'Sedgeford_Mean';

% savename
S.savename   =  'Sedgeford_Mean_tf_FootHeight_ry02';

% impose the stride frequency
S.Bounds.tf = [FinalTime-0.01 FinalTime+0.01];

% impose a minimal foot height during swing
S.Constr.ImposeFootHeight = true; % imposing contraint on foot height during swing
S.Constr.FootHeight.Vx_Treshold = S.v_tgt; % minimal forward velocity of the foot when constraint on footheight is applied
S.Constr.FootHeight.Ry_Treshold = 0.15; % minimal height of the foot (when velocity constraint is active)
S.Constr.FootHeight.b = 10;
% run the optimization
f_PredSim_Gait92(S);     % run the optimization
f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-proces simulation results
