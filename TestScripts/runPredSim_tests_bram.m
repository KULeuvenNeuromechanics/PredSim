
% current idea
[S] = getDefaultSettings(S)

[model_info,S] = preProcess(S,model_info)

    [~]          = osim2dll(osim_model,S)
    
    [model_info] = getModelInfo(osim_model, S)
    
    [model_info] = scaleMTParams(S,model_info)
    
    [~] = MuscleAnalysis(model_info,S)
        
        [] = getBounds
        
            [] = getIK