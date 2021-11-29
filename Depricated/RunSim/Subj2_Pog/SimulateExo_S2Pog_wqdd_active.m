%% Default settings

clear all; close all; clc;
% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 8;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% output folder
S.ResultsFolder     = 'S2_SensQdd';

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

S.Bounds.ActLower = 0.01;

S.W.E       = 500;      % weight metabolic energy rate
S.W.Ak      = 50000;    % weight joint accelerations
S.W.ArmE    = 10^6;     % weight arm excitations
S.W.passMom = 1000;     % weight passive torques
S.W.A       = 2000;     % weight muscle activations
S.W.exp_E   = 2;        % power metabolic energy
S.W.Mtp     = 10^6;     % weight mtp excitations
S.W.u       = 0.001;    % weight on excitations arms actuators
S.W.Lumbar  = 10^5;
    
    

% Simulation without exoskeleton
QddVect = [1000 10000 20000 50000 80000 100000 200000];

for i=1:length(QddVect)    
    S.W.Ak          = QddVect(i);    % weight joint accelerations
    S.ExoBool       = 1;    
    S.DataSet       = 'Poggensee_2020_Subj2_T1';
    S.ExoScale      = 1;
    S.ExternalFunc  = 'SimExo_Pog_s2_talus_out.dll';        % external function
    S.ExternalFunc2 = 'SimExo_Pog_s2_ExportAll.dll';        % external function for post processing
    S.savename      = ['Active_T1_qdd' num2str(S.W.Ak)];
    f_PredSim_Gait92(S);     % run the optimization    
end



