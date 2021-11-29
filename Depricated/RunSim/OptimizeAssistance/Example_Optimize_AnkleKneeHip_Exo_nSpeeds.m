%% Example optimize ankle-foot exoskeleton assistance

clear all; close all; clc;
% settings for optimization
S.N         = 50;       % number of mesh intervals
S.NThreads  = 8;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% Folder with default functions
S.subject            = 's1_Poggensee';

% output folder
S.ResultsFolder     = 'Example_ExoDesign_Speeds';

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
S.OptTexo_AnkleKneeHip.Bool = true;
S.OptTexo_AnkleKneeHip.Tbound_Ankle = [-100 0]; 
S.OptTexo_AnkleKneeHip.Tbound_Knee = [-100 100]; 
S.OptTexo_AnkleKneeHip.Tbound_Hip = [-100 100]; 
S.OptTexo_AnkleKneeHip.Pbound_Ankle = [-3000 3000]; 
S.OptTexo_AnkleKneeHip.Pbound_Knee = [-3000 3000]; 
S.OptTexo_AnkleKneeHip.Pbound_Hip = [-3000 3000]; 


vSpeeds = [0.6:0.2:1.6];
for i=1:length(vSpeeds)
    S.v_tgt     = vSpeeds(i);     % average speed
    % Simulation with active exoskeleton
    S.ExoBool       = 1;
    S.ExoScale      = 1;
    S.savename      = ['Active_' num2str(round(vSpeeds*10))];
    S.ExternalFunc  = 'SimAnkleKneeHipExo.dll';        % external function
    S.ExternalFunc2 = 'SimAnkleKneeHipExo_pp.dll';        % external function for post processing
    f_PredSim_Gait92(S);     % run the optimization
    %     f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-proces simulation results
end
for i=1:length(vSpeeds)
    S.v_tgt     = vSpeeds(i);     % average speed
    % Simulation with passive exoskeleton
    S.ExoBool       = 0;
    S.ExoScale      = 0;
    S.OptTexo_AnkleKneeHip.Bool = false;
    S.savename      = ['Passive_' num2str(round(vSpeeds*10))];
    S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % external function
    S.ExternalFunc2 = 'SimExo_3D_ExportAll.dll';        % external function for post processing
    f_PredSim_Gait92(S);     % run the optimization
    %     f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-proces simulation results
end
