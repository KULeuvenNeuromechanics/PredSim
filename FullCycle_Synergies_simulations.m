max_iter = 4000;

%% 5 + 5
Syn = 1;
NSyn_r = 5;
NSyn_l = 5;
weights.Syn_constr = 1e4; % cost function weight for (a-WH)^2
SynConstrLower = -0.001;
SynConstrUpper = 0.001;
gaitmotion_type = 'HalfGaitCycle';
% gaitmotion_type = 'FullGaitCycle';
sim_name = 'Syn_5R_5L_HalfCycle_onlyMeshPoints';

main;

if(false)
    
    %% 4 + 4
    Syn = 1;
    NSyn_r = 4;
    NSyn_l = 4;
    weights.Syn_constr = 1e4; % cost function weight for (a-WH)^2
    SynConstrLower = -0.001;
    SynConstrUpper = 0.001;
    gaitmotion_type = 'FullGaitCycle';
    sim_name = 'Syn_4R_4L_FullCycle';
    
    main;
    
    %% 3 + 3
    Syn = 1;
    NSyn_r = 3;
    NSyn_l = 3;
    weights.Syn_constr = 1e4; % cost function weight for (a-WH)^2
    SynConstrLower = -0.001;
    SynConstrUpper = 0.001;
    gaitmotion_type = 'FullGaitCycle';
    sim_name = 'Syn_3R_3L_FullCycle';
    
    main;
    
    %% 2 + 2
    Syn = 1;
    NSyn_r = 2;
    NSyn_l = 2;
    weights.Syn_constr = 1e4; % cost function weight for (a-WH)^2
    SynConstrLower = -0.001;
    SynConstrUpper = 0.001;
    gaitmotion_type = 'FullGaitCycle';
    sim_name = 'Syn_2R_2L_FullCycle';
    
    main;
    
    %% 5 + 4
    Syn = 1;
    NSyn_r = 5;
    NSyn_l = 4;
    weights.Syn_constr = 1e4; % cost function weight for (a-WH)^2
    SynConstrLower = -0.001;
    SynConstrUpper = 0.001;
    gaitmotion_type = 'FullGaitCycle';
    sim_name = 'Syn_5R_4L_FullCycle';
    
    main;
    
    %% 5 + 3
    Syn = 1;
    NSyn_r = 5;
    NSyn_l = 3;
    weights.Syn_constr = 1e4; % cost function weight for (a-WH)^2
    SynConstrLower = -0.001;
    SynConstrUpper = 0.001;
    gaitmotion_type = 'FullGaitCycle';
    sim_name = 'Syn_5R_3L_FullCycle';
    
    main;
    
    %% 5 + 2
    Syn = 1;
    NSyn_r = 5;
    NSyn_l = 2;
    weights.Syn_constr = 1e4; % cost function weight for (a-WH)^2
    SynConstrLower = -0.001;
    SynConstrUpper = 0.001;
    gaitmotion_type = 'FullGaitCycle';
    sim_name = 'Syn_5R_2L_FullCycle';
    
    main;
    
end