function [] = run_pred_sim(S,osim_path)

addpath([S.misc.main_path '\VariousFunctions'])

% Need to add casadi to path inside this function when running it as batch
% job. Passing the folder as AdditionalPaths when creating the job does not
% work. Feel free to change if you find a cleaner solution.
if S.solver.run_as_batch_job
    addpath(genpath(S.solver.CasADi_path));
end

S = getDefaultSettings(S);


%% PreProcessing
addpath([S.misc.main_path '\PreProcessing'])
disp('Start PreProcessing...')
t0 = tic;
[S,model_info] = PreProcessing(S,osim_path);
disp(['... PreProcessing done. Time elapsed ' num2str(toc(t0)) ' s'])

%% Creating casadi functions
addpath([S.misc.main_path '\CasadiFunctions'])
disp('Start creating CasADi functions...')
t0 = tic;
[f_casadi] = createCasadiFunctions(S,model_info);
disp(['... CasADi functions created. Time elapsed ' num2str(toc(t0)) ' s'])

%% Formulating OCP
addpath([S.misc.main_path '\OCP'])
OCP_formulation(S,model_info,f_casadi);

%% PostProcessing
addpath([S.misc.main_path '\PostProcessing'])
PostProcessing(S,model_info,f_casadi);
