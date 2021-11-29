function [S] = GetDefaultSettings(S)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%% Default settings

if ~isfield(S,'linear_solver')    
    S.linear_solver = 'mumps';
end

if ~isfield(S,'tol_ipopt')
    S.tol_ipopt     = 4;
end

if ~isfield(S,'savename_ig')
    S.savename_ig= [];
end

if ~isfield(S,'ResultsF_ig')
    S.ResultsF_ig = [];
end

if  ~isfield(S,'NThreads') || isempty(S.NThreads)
    S.NThreads = 4;
end

if ~isfield(S,'mass') || isempty(S.mass)
    S.mass = 64;
end

if ~isfield(S,'subject') || isempty(S.subject)
    S.subject = 'subject1';
end

% quasi random initial guess
if ~isfield(S,'IG_PelvisY') || isempty(S.IG_PelvisY)
    S.IG_PelvisY = 0.9385;
end

% default settings walking speed
if ~isfield(S,'v_tgt') || isempty(S.v_tgt)
    S.v_tgt = 1.25;
end

% default number of mesh intervals
if ~isfield(S,'N') || isempty(S.N)
    S.N         = 50;       
end

% default weights
if isfield(S,'W')
    if ~isfield(S.W,'E')
        S.W.E       = 500;      % weight metabolic energy rate
    end
    if ~isfield(S.W,'Ak')
        S.W.Ak      = 50000;    % weight joint accelerations
    end
    if ~isfield(S.W,'ArmE')
        S.W.ArmE    = 10^6;     % weight arm excitations
    end
    if ~isfield(S.W,'passMom')
        S.W.passMom = 1000;     % weight passive torques
    end
    if ~isfield(S.W,'A')
        S.W.A       = 2000;     % weight muscle activations
    end
    if ~isfield(S.W,'exp_E')
        S.W.exp_E   = 2;        % power metabolic energy
    end
    if ~isfield(S.W,'Mtp')
        S.W.Mtp     = 10^6;     % weight mtp excitations
    end
    if ~isfield(S.W,'u')
        S.W.u       = 0.001;    % weight on excitations arms actuators
    end
    if ~isfield(S.W,'Lumbar')
        S.W.Lumbar  = 10^5;
    end
    if ~isfield(S.W,'EccF')
        S.W.EccF =0;
    end
else
    S.W.E       = 500;      % weight metabolic energy rate
    S.W.Ak      = 50000;    % weight joint accelerations
    S.W.ArmE    = 10^6;     % weight arm excitations
    S.W.passMom = 1000;     % weight passive torques
    S.W.A       = 2000;     % weight muscle activations
    S.W.exp_E   = 2;        % power metabolic energy
    S.W.Mtp     = 10^6;     % weight mtp excitations
    S.W.u       = 0.001;    % weight on excitations arms actuators
    S.W.Lumbar  = 10^5;     % minimize exctitation lumbar joint
    S.W.EccF    = 0;        % minimize eccentric muscle force (squared)
end

% initial guess identifier (1: quasi random, 2: data-based)
if ~isfield(S,'IGsel')
    S.IGsel     = 1;        
end

% initial guess mode identifier (1 walk, 2 run, 3prev.solution)
if ~isfield(S,'IGmodeID')
    S.IGmodeID  = 1;        
end

% initial guess case identifier
if ~isfield(S,'IGcase')    
    S.IGcase    = 0;        
end

% weakness hip actuators
if ~isfield(S,'h_weak')
    S.h_weak    = 0;     
end

% maximal contraction velocity identifier
if ~isfield(S,'Max_s')
    S.Max_s     = 0;    
end

% weakness ankle plantaflexors
if ~isfield(S,'pf_weak')
    S.pf_weak   = 0;      
end

% metabolic energy model identifier
if ~isfield(S,'mE')
    S.mE        = 0;       
end

% co-contraction identifier
if ~isfield(S,'coCont')
    S.coCont    = 0;        
end

% Kinematics Constraints - Default Settings
if isfield(S,'Constr')
    if ~isfield(S.Constr,'calcn')
        S.Constr.calcn = 0.09;  % by default at least 9cm distance between calcn
    end
    if ~isfield(S.Constr,'toes')
        S.Constr.toes = 0.1; % by default at least 10cm distance between toes
    end
    if ~isfield(S.Constr,'tibia')
        S.Constr.tibia = 0.11; % by default at least 11cm distance between toes
    end
    if ~isfield(S.Constr,'ImposeFootHeight')
        S.Constr.ImposeFootHeight = false;
    end
    if ~isfield(S.Constr,'FootHeight') || ~isfield(S.Constr.FootHeight,'Vx_Treshold')
        S.Constr.FootHeight.Vx_Treshold = S.v_tgt;
    end
    if ~isfield(S.Constr,'FootHeight') || ~isfield(S.Constr.FootHeight,'Ry_Treshold')
        S.Constr.FootHeight.Ry_Treshold = 0.1;
    end
    if ~isfield(S.Constr,'FootHeight') || ~isfield(S.Constr.FootHeight,'b')
        S.Constr.FootHeight.b = 10;
    end
else
    S.Constr.calcn = 0.09;  % by default at least 9cm distance between calcn
    S.Constr.toes = 0.1; % by default at least 10cm distance between toes
    S.Constr.tibia = 0.11; % by default at least 11cm distance between toes
    S.Constr.ImposeFootHeight = false; % imposing contraint on foot height during swing
    S.Constr.FootHeight.Vx_Treshold = S.v_tgt; % minimal forward velocity of the foot when constraint on footheight is applied
    S.Constr.FootHeight.Ry_Treshold = 0.1; % minimal height of the foot (when velocity constraint is active)
    S.Constr.FootHeight.b = 10; % smoothing factor tanh foot height constraint (velocity part)
end

% Settings related to bounds on muscle activations
if isfield(S,'Bounds')
    if ~isfield(S.Bounds,'ActLower')
        S.Bounds.ActLower = 0.02;
    end
    if ~isfield(S.Bounds,'ActLowerHip')
        S.Bounds.ActLowerHip = [];
    end
    if ~isfield(S.Bounds,'ActLowerKnee')
        S.Bounds.ActLowerKnee = [];
    end
    if ~isfield(S.Bounds,'ActLowerAnkle')
        S.Bounds.ActLowerAnkle = [];
    end    
else
    S.Bounds.ActLower = 0.02;
    S.Bounds.ActLowerHip = [];
    S.Bounds.ActLowerKnee = [];
    S.Bounds.ActLowerAnkle = [];
end

% different function for bounds
if ~isfield(S,'Bounds_Running')
    S.Bounds_Running = false;
end

% bounds on final time (i.e. imposing stride frequency / stride length)
if ~isfield(S.Bounds,'tf')
    S.Bounds.tf = [];
end

% scaling exoskeleton timing
if isfield(S,'PercStance')
    if ~isfield(S.PercStance,'bool')
        S.PercStance.bool = 0;        
    end
    if ~isfield(S.PercStance,'xStanceOr')
        S.PercStance.xStanceOr = 0.61; % duration stance phase in experiment
    end
    if ~isfield(S.PercStance,'xStanceNew')
        S.PercStance.xStanceNew = 0.58; % duration stance phase in simulation
    end
else
    S.PercStance.bool = 0;
end

% parallel computation settings
if ~isfield(S,'parallelMode')
    S.parallelMode = 'thread';
end

% symmetric motion ?
if ~isfield(S,'Symmetric')
    S.Symmetric = true;
end

% oeriodic motion
if ~isfield(S,'Periodic')
    S.Periodic = false;
end

% model selection
if ~isfield(S,'ModelName')
    S.ModelName = 'Gait92';
end

% default IK file to determine bounds
if ~isfield(S,'IKfile_Bounds')
    S.IKfile_Bounds = 'OpenSimModel\IK_Bounds_Default.mat';
end

% default IK file for initial guess (when used data-informed guess)
if ~isfield(S,'S.IKfile_guess')
    S.IKfile_guess = 'OpenSimModel\IK_Guess_Default.mat';
end

% path with exoskeleton torque profile
if ~isfield(S,'DataSet')
    S.DataSet = 'PoggenSee2020_AFO';
end

% Boolean for exoskeleton use
if ~isfield(S,'ExoBool')
    S.ExoBool       = 0;
end

% scaling assistance profile
if ~isfield(S,'ExoScale')
    S.ExoScale      = 0;        % scale factor of exoskeleton assistance profile = 0 (i.e. no assistance)
end

% design assistance profile - Ankle
if ~isfield(S,'OptTexo_Ankle')
    S.OptTexo_Ankle.Bool = false;
end

% design assistance profile - Knee
if ~isfield(S,'OptTexo_Kee')
    S.OptTexo_Knee.Bool = false;
end

% design assistance profile - Hip
if ~isfield(S,'OptTexo_Hip')
    S.OptTexo_Hip.Bool = false;
end

% design assistance profile - Ankle, Knee, Hip
if ~isfield(S,'OptTexo_AnkleKneeHip')
    S.OptTexo_AnkleKneeHip.Bool = false;
elseif  S.OptTexo_AnkleKneeHip.Bool
    if ~isfield(S.OptTexo_AnkleKneeHip,'Tbound_Ankle')
        S.OptTexo_AnkleKneeHip.Tbound_Ankle = [-50 50];
    end
    if ~isfield(S.OptTexo_AnkleKneeHip,'Tbound_Knee')
        S.OptTexo_AnkleKneeHip.Tbound_Knee = [-50 50];
    end
    if ~isfield(S.OptTexo_AnkleKneeHip,'Tbound_Hip')
        S.OptTexo_AnkleKneeHip.Tbound_Hip = [-50 50];
    end
end


%% objective function
if ~isfield(S,'EModel')
    % options are
    % (1) Bhargave2014
    % (2) Marg1968
    S.EModel = 'Bhargava2014';
end



%% Booleans for flow control based on input
if S.OptTexo_Ankle.Bool || S.OptTexo_Knee.Bool || S.OptTexo_Hip.Bool || ...
        S.OptTexo_AnkleKneeHip.Bool
    S.ExoDesignBool = true;
else
    S.ExoDesignBool = false;
end



% Print the settings to the screen
disp(S);





end

