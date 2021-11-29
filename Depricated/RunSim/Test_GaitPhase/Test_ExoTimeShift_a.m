%% Batch Run StepWidth

% Influence of the kinematic constraint on stepwidth
clear all; close all; clc;

%% Default settings

% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 4;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% Folder with default functions
S.subject            = 's1_Poggensee';

% output folder
S.ResultsFolder     = 'Sens_TimeShift';

% select tendon stiffness of 20
S.CasadiFunc_Folders = 'Casadi_s1Pog_ScaleParam_k20';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)
S.savename_ig   = 'NoExo';

%% Imposing change in stepwidth with passive exoskeleton

S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % this one is with the pinjoint mtp
S.ExoBool       = 1;
S.ExoScale      = 1;
S.DataSet       = 'PoggenSee2020_AFO';
S.ExoTiming.Stance = 56;
S.savename      = 'TimingOffset_Stance56';

% to the Run Sim path
CurrentPath = pwd;
[RunSimPath,~] = fileparts(CurrentPath);
cd(RunSimPath);

% Run with the time shift
% f_PredSim_Pog_CalcnT_ShiftExoSupport(S);

% sensitivity
for i = 55:1:59
    S.ExoTiming.Stance = i;
    S.savename      = ['TimingOffset_Stance' num2str(i)];
    f_PredSim_Pog_CalcnT_ShiftExoSupport(S);
end


% Run without the time shift
S.savename      = 'NoShift';
f_PredSim_PoggenSee2020_CalcnT(S);

% Back to origin al path
cd(CurrentPath);
