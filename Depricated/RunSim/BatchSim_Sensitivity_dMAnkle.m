%% Default settings

clear all; close all; clc;
% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 8;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% Folder with default functions
S.subject            = 's1_Poggensee';

% output folder
S.ResultsFolder     = 'Sensitivity_dMAnkle';

% select folder with polynomials
S.PolyFolder = 's1_Poggensee';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3 prev.solution, 4 solution from /IG/Data folder)
S.savename_ig   = 'NoExo';  % name of the IG (.mot) file

% dataset with exoskeleton torque profile
S.DataSet       = 'PoggenSee2020_AFO';

% Select model
S.ModelName = 'Gait92'; % other option is 'Rajagopal'

% mass of the subject
S.mass = 64;    % mass in kg

%% Loop over folders with Ankle dM
% (casadifunctsion generated with Sensiity_dMAnkle_Gait92.m)

dM_Vect = -0.02:0.005:0.02;
for i=1:length(dM_Vect)
    % select moment arm offset
    dMSel = dM_Vect(i);
    
    % Folder with CasadiFunctions    
    dM_string = num2str(round(dMSel*1000));
    S.CasadiFunc_Folders = ['s1_Poggensee_dMAnkle_' dM_string];  

    % Simulation without exoskeelton
    S.ExoBool       = 0;
    S.ExoScale      = 0;        % scale factor of exoskeleton assistance profile = 0 (i.e. no assistance)
    S.ExternalFunc  = 'PredSim_3D_Pog_s1_mtp.dll';        % external function
    S.ExternalFunc2 = 'PredSim_3D_Pog_s1_mtp_pp.dll';     % external function for post-processing
    S.savename      = ['dMAnkle_' dM_string];
    f_PredSim_Gait92(S);     % run the optimization
    f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-proces simulation results
end
