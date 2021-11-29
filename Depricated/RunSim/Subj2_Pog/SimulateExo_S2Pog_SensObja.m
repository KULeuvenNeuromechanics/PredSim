%% Default settings

clear all; close all; clc;
% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 2;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% output folder
S.ResultsFolder     = 'S2_PogSensW';

% selection folder with Casadi Functions
S.CasadiFunc_Folders = 'S2_Pog'; 

% select folder with polynomials
S.PolyFolder = 'S2_Pog';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3 prev.solution, 4 solution from /IG/Data folder)
S.savename_ig   = 'NoExo';  % name of the IG (.mot) file

% Select model
S.ModelName = 'Gait92'; % other option is 'Rajagopal'

% mass of the subject
S.mass = 74;    % mass in kg

S.Bounds.ActLower = 0.01; % lower bound muscle activity

% change weight objective function
S.W.E       = 1000;      % weight metabolic energy rate
S.W.Ak      = 10000;    % weight joint accelerations
S.W.ArmE    = 10^6;     % weight arm excitations
S.W.passMom = 1000;     % weight passive torques
S.W.A       = 2000;     % weight muscle activations
S.W.exp_E   = 2;        % power metabolic energy
S.W.Mtp     = 10^6;     % weight mtp excitations
S.W.u       = 0.001;    % weight on excitations arms actuators
S.W.Lumbar  = 10^5;

% Simulation without exoskeleton
S.ExoBool       = 0;    
S.ExoScale      = 0;        % scale factor of exoskeleton assistance profile = 0 (i.e. no assistance)
S.ExternalFunc  = 'PredSim_3D_Pog_s2_mtp.dll';        % external function
S.ExternalFunc2 = 'PredSim_3D_Pog_s2_mtp_pp.dll';     % external function for post-processing
S.savename      = 'NoExo_w1';
f_PredSim_Gait92(S);     % run the optimization
% f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-proces simulation results

% Simulation with passive exoskeleton
S.ExoBool       = 1;    
S.ExoScale      = 0;        % scale factor of exoskeleton assistance profile = 0 (i.e. no assistance)
S.ExternalFunc  = 'SimExo_Pog_s2_talus_out.dll';        % external function
S.ExternalFunc2 = 'SimExo_Pog_s2_ExportAll.dll';        % external function for post processing
S.savename      = 'Passive_w1';
f_PredSim_Gait92(S);     % run the optimization
% f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-proces simulation results

% Simulation with active exoskeleton
S.ExoBool       = 1;    
S.DataSet       = 'Poggensee_2020_Subj2_T2';
S.ExoScale      = 1;    
S.savename      = 'Active_T2_w1';
S.ExternalFunc  = 'SimExo_Pog_s2_talus_out.dll';        % external function
S.ExternalFunc2 = 'SimExo_Pog_s2_ExportAll.dll';        % external function for post processing
f_PredSim_Gait92(S);     % run the optimization
% f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-proces simulation results

% Simulation with active exoskeleton
S.ExoBool       = 1;    
S.DataSet       = 'Poggensee_2020_Subj2_T1';
S.ExoScale      = 1;    
S.savename      = 'Active_T1_w1';
S.ExternalFunc  = 'SimExo_Pog_s2_talus_out.dll';        % external function
S.ExternalFunc2 = 'SimExo_Pog_s2_ExportAll.dll';        % external function for post processing
f_PredSim_Gait92(S);     % run the optimization
% f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-proces simulation results




