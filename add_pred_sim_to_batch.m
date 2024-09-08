function [varargout] = add_pred_sim_to_batch(S,osim_path,varargin)
% --------------------------------------------------------------------------
% add_pred_sim_to_batch
%   This functions runs the predictive simulation as a batch job. Doing so
%   allows the simulation to run in the background. It is possible to run
%   multiple simulations at the same time, and queue more to run when they
%   are done.
%   Requires the parallel computing toolbox
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%   - PredSim_fun (optional) -
%   * handle of function that is added to the batch. By default this is
%   run_pred_sim. (Advanced)
%
% OUTPUT:
%   - S (optional)-
%   * setting structure S
% 
%
% Original author: Lars D'Hondt
% Original date: 09/May/2022
% --------------------------------------------------------------------------

if isempty(varargin)
    PredSim_fun = 'run_pred_sim';
else
    PredSim_fun = varargin{1};
end
% Running a function as part of batch needs all paths that are called inside
% the function
additional_paths = {};
[pathOsim,~,~] = fileparts(osim_path);
additional_paths{end+1} = pathOsim;
additional_paths{end+1} = fullfile(S.misc.main_path,'DefaultSettings');
additional_paths{end+1} = fullfile(S.misc.main_path,'PreProcessing');
additional_paths{end+1} = fullfile(S.misc.main_path,'CasadiFunctions');
additional_paths{end+1} = fullfile(S.misc.main_path,'OCP');
additional_paths{end+1} = fullfile(S.misc.main_path,'Subjects');
additional_paths{end+1} = fullfile(S.misc.main_path,'VariousFunctions');
additional_paths{end+1} = fullfile(S.misc.main_path,'PostProcessing');

% Select parallel cluster
all_par_clusters = parallel.clusterProfiles;
if isfield(S.solver,'par_cluster_name') && ~isempty(S.solver.par_cluster_name)...
        && ismember(S.solver.par_cluster_name,all_par_clusters)
    myCluster = parcluster(S.solver.par_cluster_name);
else
    myCluster = parcluster;
end

% Adapt number of threads if needed
N_threads = myCluster.NumThreads;
if isfield(S.solver,'N_threads') && S.solver.N_threads > N_threads
    S.solver.N_threads = N_threads;
end

try myCluster.Jobs(end,1);
    S.solver.job_id = myCluster.Jobs(end,1).ID+1;
catch
    S.solver.job_id = 1;
end

% Add job to batch
batch(myCluster,PredSim_fun,0,{S,osim_path},'CurrentFolder',S.misc.main_path,...
    'AdditionalPaths',additional_paths);

% Return optional outputs
if nargout >= 1
    varargout{1} = S.solver.job_id;
end
