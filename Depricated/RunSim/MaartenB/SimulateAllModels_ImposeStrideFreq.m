clear all; close all; clc;

% settings for optimization
S.v_tgt     = 4;     % average speed
S.N         = 30;       % number of mesh intervals
S.NThreads  = 8;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% output folder
S.ResultsFolder     = 'AllModels_MaartenB_ImposeFreq';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 2;        % initial guess mode identifier (1 walk, 2 run, 3 prev.solution, 4 solution from /IG/Data folder)
S.IKfile_guess = 'OpenSimModel\IK_Guess_Default.mat'; % or OpenSimModel\IK_Guess_Default_Running.mat
S.IKfile_Bounds = 'OpenSimModel\IK_Bounds_Default.mat'; % or OpenSimModel\IK_Guess_Default_Running.mat

% Select model
S.ModelName = 'Gait92'; % other option is 'Rajagopal'

% mass of the subject
S.mass = 64;    % mass in kg (needed ?)

% Simulation without exoskeelton
S.ExternalFunc  = 'PredSim_3D_Pog_s1_mtp.dll';        % external function
S.ExternalFunc2 = 'PredSim_3D_Pog_s1_mtp_pp.dll';     % external function for post-processing

% bounds based on running
S.Bounds_Running = true;

FinalTime = [0.25 0.3 0.35 0.4 0.45];
for i=1:length(FinalTime)
    % selection folder with Casadi Functions
    S.CasadiFunc_Folders = 'Sedgeford_Mean_Fiso2';
    % select folder with polynomials
    S.PolyFolder = 'Sedgeford_Mean';
    % savename
    S.savename   =  ['Sedgeford_Mean_tf' num2str(i)];
    % impose the stride frequency
    S.Bounds.tf = [FinalTime(i)-0.01 FinalTime(i)+0.01];
    % run the optimization
    f_PredSim_Gait92(S);     % run the optimization
    f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-proces simulation results
end

