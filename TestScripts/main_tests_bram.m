clear all

[S] = initializeSettings();

load("C:\Users\u0138016\Downloads\model_info.mat");



S.bounds.a.lower = 0.02;

% S.subject.IG_selection = "quasi-random";


S.subject.save_folder = "C:\Users\u0138016\OneDrive - KU Leuven\GitHub\PredSim\Subjects\test";
S.subject.name = "Bram_test"; 
% S.subject.IG_selection = "quasi-random";
S.subject.IG_selection = "C:\Users\u0138016\OneDrive - KU Leuven\SimCP_2\Subjects\CP3\T0\IK\CP3_T0_07_IK.mot";
S.subject.IG_bounds = "C:\Users\u0138016\OneDrive - KU Leuven\SimCP_2\Subjects\CP3\T0\IK\CP3_T0_07_IK.mot";
% S.subject.IG_bounds = [];
% S.subject.mass = 45;
% S.subject.IG_pelvis_y = 1;

% S.subject.IG_selection = "quasi-random";
% S.subject.IG_selection = "C:\Users\u0138016\OneDrive - KU Leuven\SimCP_2\Subjects\CP3\T0\IK\CP3_T0_07_IK.mot";

S.subject.muscle_strength.names = ["glut_med1_r", "glut_med2_r"];
S.subject.muscle_strength.scale_factors = [0.8, 0.7]; 


model_info.mass = 50;
model_info.pelvis_y = 0.95;

runPredSim_tests_bram(S,osim_model);

[S] = getDefaultSettings(S);


% % current idea
% runPredSim(S,osim_model)
% 
%     [model_info]   = readOsim(osim_model)
% 
%     [S,model_info] = getDefaultSettings(S,model_info)
%     
%     [model_info]   = scaleMTParams(model_info)
%     
%     
% % another idea
% runPredSim(S,osim_model)
% 
%      [S]          = getDefaultSettings(S)    
%      
%      [model_info] = getModelInfo(S,osim_model)
% 
%         [model_info]   = scaleMTParams(S,default_model_info)
%     
    