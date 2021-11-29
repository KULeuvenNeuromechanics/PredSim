%% Settings
clear all; close all; clc;
% boolean to create casadifunctions
Bool_CreateCasFunc = false;
Bool_RunSim = true;

% path to 3D predictsim repository
% MainPath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\PredSim_3D';
MainPath = 'C:\Users\u0088756\Documents\Software\PredSim_3D';

%% Create the casadifunctions

% path information
% settings:
% Folder to save the polynomials
S.PolyFolder = 's1_Poggensee';
% Modelpath
S.ModelPath = fullfile(MainPath,'OpenSimModel','subject1_Poggensee_scaled.osim');
% model selection options: Rajagopal, Gait92
S.ModelName = 'Gait92';
% specific settings for exporting casadi functions
SettingsCasFunc.kTendon_CalfM = 20;
SettingsCasFunc.kMTP = 1.5/(pi/180)/5;
SettingsCasFunc.dMTP = 0.5;

% Create only the casadi functions
% create casadi functions for equations in optimiztion problem
if Bool_CreateCasFunc
    % adapt the damping value at muscle contraction level
    d_Muscle_Vect = [0.001 0.01 0.05 0.1 1];
    for i=1:length(d_Muscle_Vect)
        S.CasadiFunc_Folders = ['s1_dMuscle_' num2str(i)];
        S.dMuscle = d_Muscle_Vect(i);
        CreateCasadiFunctions(MainPath, S.ModelName, S.ModelPath, S.CasadiFunc_Folders,...
            S.PolyFolder,SettingsCasFunc);
    end

    % damping value at joint level
    d_Joint = [0.001 0.01 0.05 0.1 1];
    for i=1:length(d_Joint)
        S.CasadiFunc_Folders = ['s1_dJoint_' num2str(i)];
        S.dLimit = d_Joint(i);
        CreateCasadiFunctions(MainPath, S.ModelName, S.ModelPath, S.CasadiFunc_Folders,...
            S.PolyFolder,SettingsCasFunc);
    end

    % damping value at arms
    d_Arm = [0.01 0.05 0.1 1 10];
    for i=1:length(d_Arm)
        S.CasadiFunc_Folders = ['s1_dArm_' num2str(i)];
        Settings.dampingArm = d_Arm(i);
        CreateCasadiFunctions(MainPath, S.ModelName, S.ModelPath, S.CasadiFunc_Folders,...
            S.PolyFolder,SettingsCasFunc);
    end
end

%% Run the simulations
if Bool_RunSim

    % generic settings
    S.N         = 50;       % number of mesh intervals
    S.NThreads  = 8;        % number of threads for parallel computing

    % walking speed
    S.v_tgt     	= 1.25;

    % quasi random initial guess, pelvis y position
    S.IG_PelvisY = 0.896;   % subject 1 poggensee

    % Folder with default functions
    S.subject            = 's1_Poggensee';

    % output folder
    S.ResultsFolder     = 'Sensitivity_Damping';

    % select folder with polynomials
    S.PolyFolder = 's1_Poggensee';

    % initial guess based on simulations without exoskeletons
    S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
    S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)
    S.savename_ig   = 'NoExo';

    % dataset with exoskeleton torque profile
    S.DataSet       = 'PoggenSee2020_AFO';

    % Select model
    S.ModelName = 'Gait92'; % other option is 'Rajagopal'

    % mass of the subject
    S.mass = 68;    % mass in kg

    % set model name
    S.ModelName = 'Gait92';

    % external functions
    S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % external function
    S.ExternalFunc2 = 'SimExo_3D_ExportAll.dll';        % external function for post-processing

    d_Muscle_Vect = [0.001 0.01 0.05 0.1 1];
    for i=1:length(d_Muscle_Vect)
        S.CasadiFunc_Folders = ['s1_dMuscle_' num2str(i)];
        S.ExoBool       = 1;
        S.ExoScale      = 0;
        S.savename      = ['P_s1_dMuscle_' num2str(i)];
        f_PredSim_Gait92(S);     % run the optimization

        S.ExoBool       = 1;
        S.ExoScale      = 1;
        S.savename      = ['A_s1_dMuscle_' num2str(i)];
        f_PredSim_Gait92(S);     % run the optimization
    end

%     % damping value at joint level
%     d_Joint = [0.001 0.01 0.05 0.1 1];
%     for i=1:length(d_Joint)
%         S.CasadiFunc_Folders = ['s1_dJoint_' num2str(i)];
%         S.ExoBool       = 1;
%         S.ExoScale      = 0;
%         S.savename      = ['P_s1_dJoint_' num2str(i)];
%         f_PredSim_Gait92(S);     % run the optimization
% 
%         S.ExoBool       = 1;
%         S.ExoScale      = 1;
%         S.savename      = ['A_s1_dJoint_' num2str(i)];
%         f_PredSim_Gait92(S);     % run the optimization
%     end
% 
%     % damping value at arms
%     d_Arm = [0.01 0.05 0.1 1 10];
%     for i=1:length(d_Arm)
%         S.CasadiFunc_Folders = ['s1_dArm_' num2str(i)];
%         S.ExoBool       = 1;
%         S.ExoScale      = 0;
%         S.savename      = ['P_s1_dArm_' num2str(i)];
%         f_PredSim_Gait92(S);     % run the optimization
% 
%         S.ExoBool       = 1;
%         S.ExoScale      = 1;
%         S.savename      = ['A_s1_dArm_' num2str(i)];
%         f_PredSim_Gait92(S);     % run the optimization
%     end
end



