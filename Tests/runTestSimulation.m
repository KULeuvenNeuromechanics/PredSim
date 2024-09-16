function [] = runTestSimulation(model_name, options)
% --------------------------------------------------------------------------
% runTestSimulation
%   Run a simulation and compare the outcome to a reference.
% 
%
% INPUT:
%   - model_name -
%   * Name of OpenSim model that is to be used for testing.
% 
%   - options (optional) -
%   * Options to configure the test. [name-value pairs]
%   List of options:
%       - run_as_batch_job  see S.solver.run_as_batch_job Default is false
%       - CasADi_path       see S.solver.CasADi_path Default is []
%       - OpenSimAD_compiler  see S.OpenSimADOptions.compiler Default is []
%       - OpenSimAD_verbose_mode see S.OpenSimADOptions.verbose_mode
%       - max_iter          see S.solver.max_iter Default is 10000
%       - mode  'paper': use settings from the paper of the model, 
%               'fast': use settings for faster simulation Default is 'paper'
%       - compare_printout  compare ipopt printout. Default is true. Needs
%       only a few iterations, but can give false negatives.
%       - compare_result    compare result. Default is true. Reliable, but
%       need simulation to run untill convergence.
%
%
% OUTPUT:
%   - (no outputs) -
% 
% Original author: Lars D'Hondt
% Original date: 4/January/2024
% --------------------------------------------------------------------------

arguments
    model_name char
    options.run_as_batch_job (1,1) logical = false;
    options.CasADi_path char = [];
    options.OpenSimAD_compiler char = [];
    options.OpenSimAD_verbose_mode (1,1) logical = false;
    options.max_iter (1,1) double = 10000;
    options.mode char {mustBeMember(options.mode,{'paper','fast'})} = 'paper';
    options.compare_printout (1,1) logical = true;
    options.compare_result (1,1) logical = true;
end


%% Configure settings for simulation

[pathTestDir,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathTestDir);

cd(pathRepo);

pathDefaultSettings = fullfile(pathRepo,'DefaultSettings');
addpath(pathDefaultSettings)
addpath(fullfile(pathRepo,'VariousFunctions'))

% load reference result
res_ref = load(fullfile(pathTestDir,'ReferenceResults',model_name,[model_name,'_',options.mode,'.mat']),'R','model_info');
R_ref = res_ref.R;


% Initialize S
S = R_ref.S;


% Store result of test simulation in Debug folder
if ~isfolder(fullfile(pathRepo,'Debug'))
    mkdir(fullfile(pathRepo,'Debug'));
end
if ~isfolder(fullfile(pathRepo,'Debug','TestResults'))
    mkdir(fullfile(pathRepo,'Debug','TestResults'));
end
S.misc.save_folder  = fullfile(pathRepo,'Debug','TestResults'); 
result_filename = [model_name '_test_' datestr(datetime,30)];
S.misc.result_filename = result_filename;

% update paths
S.solver.IG_selection = replace(S.solver.IG_selection, S.misc.main_path, pathRepo);
osim_path = replace(res_ref.model_info.osim_path, S.misc.main_path, pathRepo);

S.misc.main_path = pathRepo;


% Make sure casadi path is set up correctly
if isempty(options.CasADi_path)
    try
        S.solver.CasADi_path = casadi.GlobalOptions.getCasadiPath();
    catch
        error("Please add CasADi to the matlab search path, or pass the path " + ...
            "to your CasADi installation (top folder) to 'CasADi_path'.")
    end
elseif  ~isfolder(options.CasADi_path)
    error("Unable to find the path assigned 'CasADi_path':" + ...
        " \n\t%s",options.CasADi_path)
end

if ~isempty(options.OpenSimAD_compiler)
    S.OpenSimADOptions.compiler = options.OpenSimAD_compiler;
else
    S.OpenSimADOptions = rmfield(S.OpenSimADOptions,'compiler');
end
if options.max_iter > 0
    S.solver.max_iter = options.max_iter;
end

S.OpenSimADOptions.verbose_mode = options.OpenSimAD_verbose_mode;

%% run test simulation

if options.run_as_batch_job
    jobRunTest =@(S,osim_path) runTest(S,osim_path, model_name,options.mode,result_filename,options.compare_result,options.compare_printout,pathTestDir);
    add_pred_sim_to_batch(S,osim_path,jobRunTest);
else
    runTest(S,osim_path, model_name,options.mode,result_filename,options.compare_result,options.compare_printout,pathTestDir);
end

cd(pathTestDir);

end

%%
function [] = runTest(S,osim_path, model_name,mode,result_filename,compare_result,compare_printout,pathTestDir)


% try
    run_pred_sim(S,osim_path);
% catch
%     diary off
% end


% load reference result
res_ref = load(fullfile(pathTestDir,'ReferenceResults',model_name,[model_name,'_',mode,'.mat']),'R','model_info');
R_ref = res_ref.R;
res_test = load(fullfile(S.misc.main_path,'Debug','TestResults',[result_filename,'.mat']),'R','stats');
R_test = res_test.R;

n_iter_test = res_test.stats.iter_count;

test_success = true;

%% compare simulation results
if compare_result

    if res_test.stats.success
        fprintf('Test simulation converged\n');
    
        fields = {'kinematics.Qs_rad','kinematics.Qdots_rad','kinematics.Qddots_rad',...
            'muscles.a','muscles.da','muscles.FTtilde','muscles.dFTtilde',...
            'torque_actuators.a','torque_actuators.e','objective.absoluteValues'};
        
        max_diff = nan(length(fields),1);
        for i=1:length(fields)
            idx_dot = strfind(fields{i},'.');
            fieldname1 = fields{i}(1:idx_dot-1);
            fieldname2 = fields{i}(idx_dot+1:end);
            max_diff(i) = max(abs( R_test.(fieldname1).(fieldname2) - R_ref.(fieldname1).(fieldname2) ), [], 'all');
        end
    
        idx = find(max_diff > 10^(-S.solver.tol_ipopt));
        if isempty(idx)
            fprintf('Results of test simulation match reference up to tolerance (max difference = %f)\n',max(max_diff));
        else
            test_success = false;
            fprintf('Difference between test results and reference is above tolerance\n');
            for i=idx'
                fprintf('\t%s\tmax difference = %f\n',fields{i},max(max_diff(i)));
            end
        end

    else
        fprintf('Test simulation did not converge (%s)\n',res_test.stats.return_status);
    end
end

%% compare solver printout (logfile)
if compare_printout
    log_ref = fopen(fullfile(pathTestDir,'ReferenceResults',model_name,[model_name,'_',mode,'_log.txt']),'r');
    log_test = fopen(fullfile(S.misc.main_path,'Debug','TestResults',[result_filename,'_log.txt']),'r');
    
    try
        iter_label = 'iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls';
        
        % read lines until we encounter the iteration printout label
        isIter = false;
        while ~isIter
            next_line = fgetl(log_ref);
            isIter = strcmp(next_line, iter_label);
        end
        isIter = false;
        while ~isIter
            next_line = fgetl(log_test);
            isIter = strcmp(next_line, iter_label);
        end
        
        % loop over iterations
        testEqualsRef = true;
        nextLine_ref = fgetl(log_ref);
        nextLine_test = fgetl(log_test);
        for i=0:n_iter_test

            dat_nextLine_ref = sscanf(nextLine_ref,'%i%*c %f %f %f %f %f %*s %f %f%*c %i',9);
            dat_nextLine_test = sscanf(nextLine_test,'%i%*c %f %f %f %f %f %*s %f %f%*c %i',9);

            testEqualsRef = (dat_nextLine_ref == dat_nextLine_test);

            isNextIter = false;
            % loop over multiple lines printed within one iteration
            while ~isNextIter
                
                nextLine_ref = fgetl(log_ref);
                nextLine_test = fgetl(log_test);
        
                nextLineTrimmed = strtrim(nextLine_test);
                isNextIter = ~isempty(nextLineTrimmed) && isstrprop(nextLineTrimmed(1),'digit') || i==n_iter_test;
            end
        
            if ~testEqualsRef
                test_success = false;
                fprintf('IPOPT printouts do not match for iteration %i\n',i);
                break
            end
        end
        
        fclose(log_ref);
        fclose(log_test);
    
    catch err
        fclose(log_ref);
        fclose(log_test);
        rethrow(err)
    end
end

%%

if test_success
    fprintf(['\n=====================================================\n',...
        '|\t\t\t\t\tTest successful\t\t\t\t\t|',...
        '\n=====================================================\n\n']);

else
    f_geom_path = fullfile(S.misc.subject_path, S.misc.msk_geom_name);
    f_ID_path = fullfile(S.misc.subject_path, S.misc.external_function);

    error(['\nTest failed\n\tTry removing (or renaming or moving) ',...
        '"%s" and "%s" and rerun the test.\n'], f_geom_path, f_ID_path)
end


end