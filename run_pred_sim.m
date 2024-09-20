function [varargout] = run_pred_sim(S,osim_path)
% --------------------------------------------------------------------------
% run_pred_sim
%   This functions calls the subfunction for each step in the simulation
%   workflow.
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%
% OUTPUT:
%   - savename (optional) -
%   * results of the simulation will be saved in a file with this name
% 
% Original author: Lars D'Hondt
% Original date: March-October/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

addpath([S.misc.main_path '\VariousFunctions'])
addpath([S.misc.main_path '\WearableDevices'])

% Settings that are not specified get their default value
S = getDefaultSettings(S,osim_path);

% Add CasADi to the path
if ~isempty(S.solver.CasADi_path)
    addpath(genpath(S.solver.CasADi_path));
end

% Make sure folder to save results exists
OutFolder = S.misc.save_folder;
if ~isfolder(OutFolder)
    mkdir(OutFolder);
end

if S.post_process.load_prev_opti_vars
    % load settings and model_info when only running post-processing
    Outname = fullfile(S.misc.save_folder,[S.misc.result_filename '.mat']);
    load(Outname,'R','model_info');
    S = R.S;
    S.post_process.load_prev_opti_vars = 1;
    osim_path = model_info.osim_path;
    S = getDefaultSettings(S, osim_path); % to fill in any missing settings
    R.S = S;

elseif S.post_process.rerun
    % load settings and model_info when only running post-processing
    Outname = fullfile(S.misc.save_folder,[S.misc.result_filename '.mat']);
    load(Outname,'R','model_info');
    S = R.S;
    S.post_process.rerun = 1;
    osim_path = model_info.osim_path;
    S = getDefaultSettings(S, osim_path); % to fill in any missing settings
    
    R.S = S;

elseif isempty(S.misc.result_filename)
    if strcmp(S.misc.savename,'structured')
        % use a structured savename
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
        S.misc.result_filename = result_filename;

    elseif strcmp(S.misc.savename,'datetime')
        % use system date and time
        S.misc.result_filename = [S.subject.name '_' datestr(datetime,30)];
        
    end   
end

if nargout == 1
    varargout{1} = S.misc.result_filename;
end

%% Start diary
t00 = tic;
Outname = fullfile(S.misc.save_folder,[S.misc.result_filename '_log.txt']);
diary(Outname);
disp(' ')
disp(['Subject name: ' S.subject.name])
disp(['OpenSim model: ' osim_path])
disp(' ')
disp(' ')

%% PreProcessing
addpath([S.misc.main_path '\PreProcessing'])
disp('Start PreProcessing...')
disp(' ')
t0 = tic;
[S,model_info] = PreProcessing(S,osim_path);
disp(' ')
disp(['...PreProcessing done. Time elapsed ' num2str(toc(t0),'%.2f') ' s'])
disp(' ')
disp(' ')

%% Creating casadi functions
addpath([S.misc.main_path '\CasadiFunctions'])
addpath([S.misc.main_path '\ModelComponents'])

disp('Start creating CasADi functions...')
disp(' ')
t0 = tic;
[f_casadi] = createCasadiFunctions(S,model_info);
disp(' ')
disp(['...CasADi functions created. Time elapsed ' num2str(toc(t0),'%.2f') ' s'])
disp(' ')
disp(' ')

%% Formulating OCP
addpath([S.misc.main_path '\OCP'])
if ~S.post_process.rerun
    OCP_formulation(S,model_info,f_casadi);
    disp(' ')
    disp(' ')
end

%% PostProcessing
addpath([S.misc.main_path '\PostProcessing'])
disp('Start PostProcessing...')
disp(' ')
t0 = tic;
PostProcessing(S,model_info,f_casadi);
disp(' ')
disp(['...PostProcessing done. Time elapsed ' num2str(toc(t0),'%.2f') ' s'])
disp(' ')
disp(' ')

%% Conclude diary
disp(['Total time elapsed ' num2str(toc(t00),'%.2f') ' s'])
disp(' ')
disp(['Diary saved as ' Outname])
disp(' ')
disp(' ')
diary off


