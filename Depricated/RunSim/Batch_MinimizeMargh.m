%% Default settings

clear all; close all; clc;
% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 2;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% Folder with default functions
S.subject            = 's1_Poggensee';

% output folder
S.ResultsFolder     = 'Example_Marg1986';

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

%% Default simulations Marg1986
% set the energy model
S.EModel = 'Marg1968';

% Simulation without exoskeelton
S.ExoBool       = 0;
S.ExoScale      = 0;        % scale factor of exoskeleton assistance profile = 0 (i.e. no assistance)
S.ExternalFunc  = 'PredSim_3D_Pog_s1_mtp.dll';        % external function
S.ExternalFunc2 = 'PredSim_3D_Pog_s1_mtp_pp.dll';     % external function for post-processing

% loop over weights
ScaleOtherV = [2 1 0.8 0.5 0.4 0.3 0.1 0.01 0.001];
for i=1:length(ScaleOtherV)
    Scale = ScaleOtherV(i);
    % set the weights in the objective function
    S.W.E       = 500;          % weight metabolic energy rate
    S.W.Ak      = 50000*Scale;  % weight joint accelerations
    S.W.ArmE    = 10^6;         % weight arm excitations
    S.W.passMom = 1000;         % weight passive torques
    S.W.A       = 2000*Scale;   % weight muscle activations
    S.W.exp_E   = 2;            % power metabolic energy
    S.W.Mtp     = 10^6;         % weight mtp excitations
    S.W.u       = 0.001;        % weight on excitations arms actuators
    S.W.Lumbar  = 10^5;
    
    % run simulations
    S.savename      = ['Marg1968_w' num2str(i)];
    f_PredSim_Gait92(S);     % run the optimization
    f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-proces simulation results
end




