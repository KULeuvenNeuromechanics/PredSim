%% Default settings

clear all; close all; clc;
% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 2;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% output folder
S.ResultsFolder     = 'Pog_Subject2';

% selection folder with Casadi Functions
S.CasadiFunc_Folders = 'Casadi_s1Pog_mtp';

% select folder with polynomials
S.PolyFolder = 's1_Poggensee';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3 prev.solution, 4 solution from /IG/Data folder)
S.savename_ig   = 'NoExo';  % name of the IG (.mot) file

% Select model
S.ModelName = 'Gait92'; % other option is 'Rajagopal'

% mass of the subject
S.mass = 64;    % mass in kg

% Simulation with active exoskeleton
S.ExoBool       = 1;    
S.ExoScale      = 1;    
S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % external function
S.ExternalFunc2 = 'SimExo_3D_ExportAll.dll';        % external function for post processing

%% torque profile 1

S.DataSet       = 'Poggensee_2020_Subj2_T2';
S.savename      = 'Active_T2';
f_PredSim_Gait92(S);     % run the optimization
f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-proces simulation results


%% torque profile 2
S.DataSet       = 'Poggensee_2020_Subj2_T1';
S.savename      = 'Poggensee_2020_Subj2_T1';
f_PredSim_Gait92(S);     % run the optimization
f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-proces simulation results
