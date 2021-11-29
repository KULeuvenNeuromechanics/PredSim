
%% Batch Sim exoskeleton assistance

% Default simulations antoine with exoskeleton assistance.
% Note currently without inertia of the exoskeleton
clear all;

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
S.IG_PelvisY = 0.903;   % subject 1 poggensee
%S.IG_PelvisY = 0.9385;

% Folder with default functions
S.CasadiFunc_Folders = 'Casadi_s1Pog_mtp';
S.subject            = 's1_Poggensee';

% output folder
S.ResultsFolder = 'BatchSim_2020_03_13';

%% open cluster
if S.Flow.solveProblem
    myCluster = parcluster('Maarten_8Cores');
end


%% No Exoskeleton

% external function 
S.ExternalFunc = 'PredSim_3D_Pog_s1_mtp.dll';        % this one is with the pinjoint mtp 
S.ExternalFunc2 = 'PredSim_3D_Pog_s1_mtp_pp.dll';    % this one is with the pinjoint mtp 

% external function
S.savename      = 'NoExo';
S.loadname      = 'NoExo';
S.ExoBool       = 1;
S.ExoScale      = 0;
S.DataSet       = 'PoggenSee2020_AFO';
jobs(1) = batch(myCluster,'f_PredSim_PoggenSee2020_ToeConstraint',0,{S});



%% Passive exoskeleton

% external function 
S.ExternalFunc = 'SimExo_3D_Pog_s1_mtp.dll';        % this one is with the pinjoint mtp 
S.ExternalFunc2 = 'SimExo_3D_Pog_s1_mtp_pp.dll';    % this one is with the pinjoint mtp 
S.savename      = 'Passive';
S.loadname      = 'Passive';
S.ExoBool       = 1;
S.ExoScale      = 0;
S.DataSet       = 'PoggenSee2020_AFO';
jobs(2) = batch(myCluster,'f_PredSim_PoggenSee2020_ToeConstraint',0,{S});

%% 100% assistance exoskeleton

% external function 
S.ExternalFunc = 'SimExo_3D_Pog_s1_mtp.dll';        % this one is with the pinjoint mtp 
S.ExternalFunc2 = 'SimExo_3D_Pog_s1_mtp_pp.dll';    % this one is with the pinjoint mtp 
S.savename      = 'Active';
S.loadname      = 'Active';
S.ExoBool       = 1;
S.ExoScale      = 1;
S.DataSet       = 'PoggenSee2020_AFO';
jobs(3) = batch(myCluster,'f_PredSim_PoggenSee2020_ToeConstraint',0,{S});

%% Loop over level of assistance
% external function 
S.ExternalFunc = 'SimExo_3D_Pog_s1_mtp.dll';        % this one is with the pinjoint mtp 
S.ExternalFunc2 = 'SimExo_3D_Pog_s1_mtp_pp.dll';    % this one is with the pinjoint mtp 
AVect = 100:20:300;
for j = 1:length(AVect)
    S.savename      = ['Active_' num2str(AVect(j))];
    S.loadname      = ['Active_' num2str(AVect(j))];
    S.ExoBool       = 1;
    S.ExoScale      = AVect(j)./100;
    S.DataSet       = 'PoggenSee2020_AFO';
    jobs(3+j) = batch(myCluster,'f_PredSim_PoggenSee2020_ToeConstraint',0,{S});
end
