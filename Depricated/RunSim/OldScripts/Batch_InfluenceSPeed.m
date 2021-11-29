
%% Batch Sim exoskeleton assistance

% Default simulations antoine with exoskeleton assistance.
% Note currently without inertia of the exoskeleton
clear all;

%% Default settings

% flow control
S.Flow.solveProblem     = 1;   % set to 1 to solve problem
S.Flow.analyseResults   = 1;   % set to 1 to analyze results
S.Flow.loadResults      = 0;   % set to 1 to load results
S.Flow.saveResults      = 1;   % set to 1 to save sens. results
S.Flow.checkBoundsIG    = 0;   % set to 1 to visualize guess-bounds
S.Flow.writeIKmotion    = 1;   % set to 1 to write .mot file

% settings for optimization
S.N         = 50;       % number of mesh intervals
S.W.E       = 500;      % weight metabolic energy rate
S.W.Ak      = 50000;    % weight joint accelerations
S.W.ArmE    = 10^6;     % weight arm excitations
S.W.passMom = 1000;     % weight passive torques
S.W.A       = 2000;     % weight muscle activations
S.W.exp_E   = 2;        % power metabolic energy
S.W.Mtp     = 10^6;     % weight mtp excitations
S.W.u       = 0.001;    % weight on excitations arms actuators
S.IGsel     = 2;        % initial guess identifier (1: quasi random, 2: data-based)
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
% S.IG_PelvisY = 0.903;   % subject 1 poggensee
% S.IG_PelvisY = 0.9385;

% external function 
S.ExternalFunc = 'PredSim_mtpPin_cm3.dll';        % this one is with the pinjoint mtp 
S.ExternalFunc2 = 'PredSim_mtpPin_pp_cm3.dll';    % this one is with the pinjoint mtp 

% Folder with default functions
S.CasadiFunc_Folders = 'Casadi_subject1_mtp2';
S.subject            = 'subject1';

%% Specific settings to runs simulation


walkspeed = 80:20:140;
for i=1:length(walkspeed)
    S.v_tgt     = walkspeed(i)./100;     % average speed
    S.savename      = ['Walk_' num2str(walkspeed(i))];
    S.loadname      = ['Walk_' num2str(walkspeed(i))];
    S.ResultsFolder = 'TestHinge_Speed';
    S.ExoBool       = 0;
    S.ExoScale      = 0;
    S.DataSet       = 'PoggenSee2020_AFO';
%     f_PredSim_PoggenSee2020_ToeConstraint(S);
%     f_LoadSim_PoggenSee2020_DefaultS(S.ResultsFolder,S.loadname);
end

%%
ColV = copper(length(walkspeed));
dPath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\Results\TestHinge_Speed';
for i = 1:length(walkspeed)
    speed = walkspeed(i);
    xLabs = 'walking speed';
    if i ==1        
        PlotResults_3DSim(fullfile(dPath,['Walk_' num2str(walkspeed(i)) '_pp.mat']),ColV(i,:),[],speed,xLabs);
    else
        PlotResults_3DSim(fullfile(dPath,['Walk_' num2str(walkspeed(i)) '_pp.mat']),ColV(i,:),gcf,speed,xLabs);
    end
end



