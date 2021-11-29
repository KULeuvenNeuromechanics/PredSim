
%% Batch Sim exoskeleton assistance

% Default simulations antoine with exoskeleton assistance.
% Note currently without inertia of the exoskeleton
clear all; close all; clc;

% add casadi to the matlab path (needed on simulation computer because
% default version is 3.4.5
rmpath(genpath('G:\PhD\Matlabcodes\Casadi\casadi'));
addpath(genpath('D:\MaartenAfschrift\softInstall\casadi'));
addpath(genpath('D:\MaartenAfschrift\GitProjects\3dpredictsim'));


%% Default settings

% flow control
S.Flow.solveProblem     = 1;   % set to 1 to solve problem
S.Flow.analyseResults   = 1;   % set to 1 to analyze results
S.Flow.loadResults      = 0;   % set to 1 to load results
S.Flow.saveResults      = 1;   % set to 1 to save sens. results
S.Flow.checkBoundsIG    = 0;   % set to 1 to visualize guess-bounds
S.Flow.writeIKmotion    = 1;   % set to 1 to write .mot file

%J_sum = 500+ 50000 + 10^6 + 1000 + 2000 + 10^6 + 0.001;

% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.W.E       = 500;      % weight metabolic energy rate
S.W.Ak      = 50000;    % weight joint accelerations
S.W.ArmE    = 10^6;     % weight arm excitations
S.W.passMom = 1000;     % weight passive torques
S.W.A       = 2000;     % weight muscle activations
S.W.exp_E   = 2;        % power metabolic energy
S.W.Mtp     = 10^6;     % weight mtp excitations
S.W.u       = 0.001;    % weight on excitations arms actuators
S.IGsel     = 1;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID  = 1;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution)
S.IGcase    = 0;        % initial guess case identifier
S.h_weak    = 0;        % weakness hip actuators
S.Max_s     = 0;        % maximal contraction velocity identifier
S.pf_weak   = 0;        % weakness ankle plantaflexors
S.mE        = 0;        % metabolic energy model identifier
S.coCont    = 0;        % co-contraction identifier
S.NThreads  = 4;        % number of threads for parallel computing

% ipopt options
S.linear_solver = 'mumps';
S.tol_ipopt     = 4;

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% Folder with default functions
S.subject            = 's1_Poggensee';

% output folder
S.ResultsFolder     = 'Batch_TendonScale_k25';

% select tendon stiffness of 20
S.CasadiFunc_Folders = 'Casadi_s1Pog_ScaleParam_k25';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 3;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution)
S.ResultsF_ig   = 'BatchSim_2020_03_17_UpdIG';  % copied from 17_03_UpdIG
S.savename_ig   = 'NoExo';

%% Passive exoskeleton




% NoExo
S.ExternalFunc  = 'PredSim_3D_Pog_s1_mtp.dll';        % this one is with the pinjoint mtp
S.ExternalFunc2 = 'PredSim_3D_Pog_s1_mtp_pp.dll';    % this one is with the pinjoint mtp
S.savename      = 'NoExo';
S.ExoBool       = 0;
S.ExoScale      = 0;
S.DataSet       = 'PoggenSee2020_AFO';
f_PredSim_PoggenSee2020_ToeConstraint(S);

% Passive exoskeletmon
S.ExternalFunc  = 'SimExo_3D_Pog_s1_mtp.dll';        % this one is with the pinjoint mtp
S.ExternalFunc2 = 'SimExo_3D_Pog_s1_mtp_pp.dll';    % this one is with the pinjoint mtp
S.savename      = 'Passive';
S.ExoBool       = 1;
S.ExoScale      = 0;
S.DataSet       = 'PoggenSee2020_AFO';
f_PredSim_PoggenSee2020_ToeConstraint(S);

% Active exoskeletmon
S.ExternalFunc  = 'SimExo_3D_Pog_s1_mtp.dll';        % this one is with the pinjoint mtp
S.ExternalFunc2 = 'SimExo_3D_Pog_s1_mtp_pp.dll';    % this one is with the pinjoint mtp
S.savename      = 'Active_100';
S.ExoBool       = 1;
S.ExoScale      = 1;
S.DataSet       = 'PoggenSee2020_AFO';
f_PredSim_PoggenSee2020_ToeConstraint(S);

% Active exoskeletmon
S.ExternalFunc  = 'SimExo_3D_Pog_s1_mtp.dll';        % this one is with the pinjoint mtp
S.ExternalFunc2 = 'SimExo_3D_Pog_s1_mtp_pp.dll';    % this one is with the pinjoint mtp
S.savename      = 'Active_160';
S.ExoBool       = 1;
S.ExoScale      = 1.6;
S.DataSet       = 'PoggenSee2020_AFO';
f_PredSim_PoggenSee2020_ToeConstraint(S);

% Active exoskeletmon
S.ExternalFunc  = 'SimExo_3D_Pog_s1_mtp.dll';        % this one is with the pinjoint mtp
S.ExternalFunc2 = 'SimExo_3D_Pog_s1_mtp_pp.dll';    % this one is with the pinjoint mtp
S.savename      = 'Active_180';
S.ExoBool       = 1;
S.ExoScale      = 1.8;
S.DataSet       = 'PoggenSee2020_AFO';
f_PredSim_PoggenSee2020_ToeConstraint(S);

% Active exoskeletmon
S.ExternalFunc  = 'SimExo_3D_Pog_s1_mtp.dll';        % this one is with the pinjoint mtp
S.ExternalFunc2 = 'SimExo_3D_Pog_s1_mtp_pp.dll';    % this one is with the pinjoint mtp
S.savename      = 'Active_210';
S.ExoBool       = 1;
S.ExoScale      = 2.1;
S.DataSet       = 'PoggenSee2020_AFO';
f_PredSim_PoggenSee2020_ToeConstraint(S);




