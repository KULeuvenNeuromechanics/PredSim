function [] = run_pred_sim(S,osim_path)

addpath([S.misc.main_path '\VariousFunctions'])
S = getDefaultSettings(S);

if ~isempty(S.solver.CasADi_path)
    addpath(genpath(S.solver.CasADi_path));
end

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
if ~S.solver.PostProcess_only
    OCP_formulation(S,model_info,f_casadi);
end

%% PostProcessing
addpath([S.misc.main_path '\PostProcessing'])
PostProcessing(S,model_info,f_casadi);
