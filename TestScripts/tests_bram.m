clear all

[S] = initializeSettings();

S.bounds.a.lower = 0.02;

S.IG.selection = "quasi-random";

S.IG.selection = "C:\Users\u0138016\OneDrive - KU Leuven\SimCP_2\Subjects\CP3\T0\IK\CP3_T0_07_IK.mot";

S.subject.muscle_strength = []


% current idea
runPredSim(S,osim_model)

    [model_info]   = readOsim(osim_model)

    [S,model_info] = getDefaultSettings(S,model_info)
    
    [model_info]   = scaleMTParams(model_info)
    
    
% another idea
runPredSim(S,osim_model)

     [S]          = getDefaultSettings(S)    
     
     [model_info] = getModelInfo(S,osim_model)

        [model_info]   = scaleMTParams(S,default_model_info)
    
    