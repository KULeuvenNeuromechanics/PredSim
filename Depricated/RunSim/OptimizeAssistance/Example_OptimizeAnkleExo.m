%% Example optimize ankle-foot exoskeleton assistance

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
S.ResultsFolder     = 'Example_ExoDesign';

% selection folder with Casadi Functions
S.CasadiFunc_Folders = 'Casadi_s1Pog_mtp';

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

% lower bound on muscle activity
S.Bounds.ActLower = 0.01;

% optimize exoskeleton assistance
S.OptTexo_Ankle.Bool = true;
S.OptTexo_Ankle.Tbound = [-50 10];      % Bound on torque
S.OptTexo_Ankle.Pbound = [ -300 300];   % Bound on power


% Simulation with active exoskeleton
S.ExoBool       = 1;    
S.ExoScale      = 1;
S.savename      = 'Active_N50';
S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % external function
S.ExternalFunc2 = 'SimExo_3D_ExportAll.dll';        % external function for post processing
f_PredSim_Gait92(S);     % run the optimization
f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-proces simulation results



