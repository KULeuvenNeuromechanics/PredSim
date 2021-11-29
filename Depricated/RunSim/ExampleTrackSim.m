%% Example TrackSim subject 1
%----------------------------
clear all; close all; clc;
% settings for optimization
S.N         = 25;       % number of mesh intervals
S.NThreads  = 4;        % number of threads for parallel computing

% Folder with default functions
S.subject            = 'Rajagopal2015';

% output folder
S.ResultsFolder     = 'Example_Rajagopal2015';

% select tendon stiffness of 20
S.CasadiFunc_Folders = 'Casadi_Rajagopal2015_LongKneeM';

% dataset with exoskeleton torque profile
S.DataSet       = 'PoggenSee2020_AFO';

% Simulation without exoskeelton
S.ExoBool       = 0;    
S.ExoScale      = 0;        % scale factor of exoskeleton assistance profile = 0 (i.e. no assistance)

% weight on lumber joint activations
S.W.Lumbar      = 10^5;

% save info
S.ResultsFolder = 'Example_Rajagopal2015_Poly';
S.savename      = 'Tracking_Rajagopal_LongKneeM';

% Tracking information
% -- weights
S.W.Qs        = 10;  % weight joint kinematics
S.W.GRF       = 1;  % weight ground reaction forces
S.W.GRM       = 1;  % weight ground reaction moments
S.W.ID_act    = 10;  % weight joint kinetics
S.W.a         = 1;  % weight muscle activations
S.W.u         = 0.001;

% select trial to track
S.Tracking.trials = {'gait_14'};
S.Tracking.subject = 'subject1';

% symmetric or periodic motion
S.Symmetric = false;
S.Periodic = false;

% select model
S.Model.Rajagopal = true;
S.Model.Default = false;

% external function
S.ExternalFunc  = 'TrackSim_Subject1.dll';        % external function
S.ExternalFunc2 = 'TrackSim_Subject1_pp.dll';     % external function for post-processing

% Run tracking simulation
f_TrackSim_Rajagopal(S);
