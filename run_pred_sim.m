function [] = run_pred_sim(S,osim_path)

addpath([S.misc.main_path '\VariousFunctions'])
S = getDefaultSettings(S);

if ~isempty(S.solver.CasADi_path)
    addpath(genpath(S.solver.CasADi_path));
end

%% Output folder
OutFolder = S.subject.save_folder;
if ~isfolder(OutFolder)
    mkdir(OutFolder);
end
if strcmp(S.post_process.savename,'structured')
    cond = 1;
    ct = 1;
    while cond
        result_filename = [S.subject.name '_v' num2str(ct)];
        if ~isfile(fullfile(OutFolder,[result_filename '.mat']))
            cond = 0;
        end
        ct = ct+1;
    end
end
if isempty(S.post_process.result_filename)
    S.post_process.result_filename = result_filename;
end


%% PreProcessing
addpath([S.misc.main_path '\PreProcessing'])
disp('Start PreProcessing...')
t0 = tic;
[S,model_info] = preprocessing(S,osim_path); %PreProcessing(S,osim_path);
disp(['... PreProcessing done. Time elapsed ' num2str(toc(t0)) ' s'])

%% Creating casadi functions
addpath([S.misc.main_path '\CasadiFunctions'])
disp('Start creating CasADi functions...')
t0 = tic;
[f_casadi] = createCasadiFunctions(S,model_info);
disp(['... CasADi functions created. Time elapsed ' num2str(toc(t0)) ' s'])

%% Formulating OCP
addpath([S.misc.main_path '\OCP'])
if ~S.post_process.rerun
    OCP_formulation(S,model_info,f_casadi);
end

%% PostProcessing
addpath([S.misc.main_path '\PostProcessing'])
PostProcessing(S,model_info,f_casadi);
