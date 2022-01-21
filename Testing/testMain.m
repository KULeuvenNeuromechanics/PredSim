%% Default settings

clear all; close all; clc;
S = initializeSettings();
% settings for optimization
S.subject.v_pelvis_x_trgt  = 1.3;         % average speed
S.solver.N_meshes  = 20;                  % number of mesh intervals
S.solver.N_threads = 2;                   % number of threads for parallel computing
S.misc.gaitmotion_type = 'HalfGaitCycle';
S.subject.mtp_type='active';

% quasi random initial guess, pelvis y position
% S.subject.IG_pelvis_y = 0.896;   % subject 1 poggensee
S.subject.IG_pelvis_y = 0.95;   % test_1

% Folder with default functions
% S.subject            = 's1_Poggensee';
S.subject.name       = 'test_1';
% mass of the subject
% S.subject.mass = 64;    % mass in kg
S.subject.mass = 75.1646;
S.Bounds_Running = false;
S.subject.IG_selection='quasi-random';
% output folder
S.ResultsFolder     = 'Example_DefaultSim';
S.subject.save_folder     = 'Example_DefaultSim';

% selection folder with Casadi Functions
% S.CasadiFunc_Folders = 'Casadi_s1Pog_mtp';

% select folder with polynomials
% S.PolyFolder = 's1_Poggensee';

% % initial guess based on simulations without exoskeletons
% S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
% S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3 prev.solution, 4 solution from /IG/Data folder)
% S.savename_ig   = 'NoExo';  % name of the IG (.mot) file
% 
% % dataset with exoskeleton torque profile
% S.DataSet       = 'PoggenSee2020_AFO';

% What is S.MuscMoAsmp?
% Select model
S.ModelName = 'Gait92'; % other option is 'Rajagopal'
S.EModel = 'Bhargava2004';

% Simulation without exoskeelton
S.ExoBool       = 0;    
S.ExoScale      = 0;        % scale factor of exoskeleton assistance profile = 0 (i.e. no assistance)
S.ExternalFunc  = 'F_test_1.dll';        % external function
% S.ExternalFunc  = 'PredSim_3D_Pog_s1_mtp.dll';        % external function
% S.ExternalFunc2 = 'PredSim_3D_Pog_s1_mtp_pp.dll';     % external function for post-processing
S.savename      = 'NoExo_out_constr';
S = getDefaultSettings_testing(S);
S = updateS(S);
S.weights.EccF = 0;
main_path = pwd;
load('model_info_2.mat')
IO = model_info.ExtFunIO;
IOfields = fields(model_info.ExtFunIO);
for i=1:length(IOfields)
    model_info.ExtFunIO.(IOfields{i}) = convert2double(model_info.ExtFunIO.(IOfields{i}));
end
[model_info] = GetIndexHelper(S,model_info);

NMuscle = length(model_info.muscle_info.params.names);
model_info.muscle_info.aTendon = 35*ones(NMuscle,1);
model_info.muscle_info.aTendon(model_info.muscle_info.IndexCalf) = 20; % Why different k for calf muslces?
model_info.muscle_info.shift = getShift(model_info.muscle_info.aTendon);

[f_casadi] = createCasadiFunctions(main_path,model_info,S);
%%
f_PredSim_Gait92(model_info,S,f_casadi);     % run the optimization
% f_LoadSim_Gait92(S.ResultsFolder,S.savename) % post-proces simulation results
