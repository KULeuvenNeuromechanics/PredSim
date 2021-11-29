
%% Batch Sim exoskeleton assistance

% Default simulations antoine with exoskeleton assistance.
% Note currently without inertia of the exoskeleton
clear all; close all; clc;

%% Default settings

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
S.ResultsFolder     = 'ExpTorque';

% select tendon stiffness of 20
S.CasadiFunc_Folders = 'Casadi_s1Pog_ScaleParam_k20';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 3;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution)
S.ResultsF_ig   = 'BatchSim_2020_03_17_UpdIG';  % copied from 17_03_UpdIG
S.savename_ig   = 'NoExo';


%% Passive exoskeleton

% Passive exoskeleton (with limited implemenation)
%I verified this, and the passive implemenation gives exactly the same
%solution as in the report 1 simulations
S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % this one is with the pinjoint mtp
S.savename      = 'Active_ExpTorque';
S.ExoBool       = 1;
S.ExoScale      = 1;
S.DataSet       = 'PoggenSee2020_Exp';
f_PredSim_PoggenSee2020_CalcnT(S);


