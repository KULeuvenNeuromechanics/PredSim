%
% outline of the structure of the functions called
%
% Authors: Dhruv Gupta, Lars D'Hondt
% Date: 21/Dec/2021
%


Main % user specified settings and osim path
    []=run_pred_sim(S,osim_path)
        S=getDefaultSettings(S) % adding settings that are not specified
        %%
        [S,model_info]=preprocessing(S,osim_path)
            []=osim2dll(S,osim_path) % save ['F_' osim_file_name '.dll'] and ['F_' osim_file_name '_IO.mat'] in subject-folder (.mat in 32bit int)
            model_info = get_model_info(S,osim_path) % load IO and create struct model_info (convert 32bit to double)
            [] = read_and_scale_MTparameters(S,osim_path,model_info) % save MTparameters in subject-folder
                MTparameters = getMTparameters(osim_path,muscle_names) % grab muscle_names from model_info
                MTparameters = scale_MTparameters(S,MTparameters) % scale lM0,vM0,lTs,FMopt,alpha0 (as specified in S)
            % only run if the geometry for the given model does not exist yet
            get_musculoskeletal_geometry_approximation(S,osim_path,model_info)
                muscle_data = muscle_analysis(S,osim_path,model_info) 
                [] = polynomial_fit(muscle_data) % save muscle_spanning_joint_INFO, muscle_info, muscle_data in subject-folder
            model_info = update_model_info() % https://github.com/KULeuvenNeuromechanics/PredSim/blob/workProgrammingRetreat/OCP/Copy_of_f_PredSim_Gait92.m#L34
        %%
        [f_casadi] = createCasadiFunctions(main_path,model_info) % number of inputs will changed based on pre-processing
        % https://github.com/KULeuvenNeuromechanics/PredSim/blob/dev-buurke/CasadiFunctions/createCasadiFunctions.m#L1
            % need to check the functions that are already in this folder
        %%
        []=OCP_formulation(S,model_info,f_casadi) 
        % https://github.com/KULeuvenNeuromechanics/PredSim/blob/workProgrammingRetreat/OCP/Copy_of_f_PredSim_Gait92.m
        % need to finish based on pre-processing and casadifunctions
            % split into 2 parts:
                % 1) common for formulation and post-processing -> run every time
                % 2) make collocation casadi formulation -> run only if solving
            % save results in specified folder
        %%
        []=post_processing(S,model_info,f_casadi)
            % load results from specified folder
            % split into parts:
                % 1) needs to always be done: extract results, assert cost, reconstruct gait cycle
                % 2) specific things
                % e.g.: [R] = getGRFS(R,S,modelinfo) to extract the GRFs from the external function
                % save R
        %%
plot_figures % plot info for a list of R

% next meeting: validation steps (early january)


































