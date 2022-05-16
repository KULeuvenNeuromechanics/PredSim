function [] = run_pred_sim(S,osim_path)

addpath([S.misc.main_path '\VariousFunctions'])

% Settings that are not specified get thier default value
S = getDefaultSettings(S);

% Add CasADi to the path
if ~isempty(S.solver.CasADi_path)
    addpath(genpath(S.solver.CasADi_path));
end

% Make sure folder to save results exists
OutFolder = S.subject.save_folder;
if ~isfolder(OutFolder)
    mkdir(OutFolder);
end

if S.post_process.rerun
    % load settings and model_info when only running post-processing
    Outname = fullfile(S.subject.save_folder,[S.post_process.result_filename '.mat']);
    load(Outname,'R','model_info');
    S = R.S;
    S.post_process.rerun = 1;
    S = getDefaultSettings(S); % to fill in any missing settings
%     osim_path = model_info.osim_path;

elseif isempty(S.post_process.result_filename)
    % use a structured savename
    if strcmp(S.post_process.savename,'structured')
        if S.solver.run_as_batch_job
            result_filename = [S.subject.name '_job' num2str(S.solver.job_id)];
        else
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
        S.post_process.result_filename = result_filename;
    end   
end

%% PreProcessing
addpath([S.misc.main_path '\PreProcessing'])
disp('Start PreProcessing...')
t0 = tic;
[S,model_info] = PreProcessing(S,osim_path); %PreProcessing(S,osim_path);
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
disp('Start PostProcessing...')
t0 = tic;
PostProcessing(S,model_info,f_casadi);
disp(['... PostProcessing done. Time elapsed ' num2str(toc(t0)) ' s'])



