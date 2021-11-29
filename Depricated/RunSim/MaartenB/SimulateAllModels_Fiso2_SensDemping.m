clear all; close all; clc;

% settings for optimization
S.v_tgt     = 4;     % average speed
S.N         = 30;       % number of mesh intervals
S.NThreads  = 8;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% output folder
S.ResultsFolder     = 'AllModels_MaartenB_JointD';

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
OutNames = {'Sedgeford_Mean'};

for i=1%:length(OutNames)
    S.ExternalFunc  = 's1Pog_dis1.dll';
    S.ExternalFunc2  = 's1Pog_dis1_pp.dll';
    % select folder with polynomials
    S.PolyFolder = OutNames{i};
    % savename
    S.savename   = OutNames{i};
    d_Vect = [0.01 0.05 0.01 0.05];
    for ij=1:length(d_Vect)
        S.CasadiFunc_Folders = ['Sedgeford_Mean_Fiso2_d_' num2str(ij)];
        S.ExoBool       = 1;
        S.ExoScale      = 0;
        S.savename      = ['Sedgeford_Mean_d_' num2str(ij)];
        f_PredSim_Gait92(S);     % run the optimization
    end
end
