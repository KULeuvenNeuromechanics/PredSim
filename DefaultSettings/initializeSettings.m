function [S] = initializeSettings(varargin)
% --------------------------------------------------------------------------
% initializeSettings
%   Initialize the settings struct
% 
% INPUT:
%   - reference_name - (optional input)
%   * pass the name of a reference model to initialize the settings
%   according to the relevant publication. 
% 
% OUTPUT:
%   - S -
%   * settings struct
%
% Original author: Bram Van Den Bosch
% Original date: 01/12/2021
% --------------------------------------------------------------------------

% get path to repo
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);

S = struct;

S.metabolicE   = [];
S.misc         = [];
S.post_process = [];
S.solver       = [];
S.subject      = [];
S.weights      = [];
S.orthosis     = [];
S.OpenSimADOptions  = [];
S.flow_control = [];

% bounds have an .upper and .lower field
S.bounds.activation_all_muscles = [];
S.bounds.SLL        = [];
S.bounds.SLR        = [];
S.bounds.dist_trav  = [];
S.bounds.t_final    = [];

% polynomial order has .lower and .upper field
S.misc.poly_order = [];

% empty array of orthoses
S.orthosis.settings = {};

% save the git hash
[S.misc.git.local_hash,S.misc.git.branch_name, S.misc.git.remote_hash] = get_git_hash(pathRepo);

% initiate for warning
S.subject.adapt_IG_pelvis_y = 0;

% save computername
if ispc
    S.misc.computername = getenv('COMPUTERNAME');
elseif isunix
    S.misc.computername = getenv('HOSTNAME');
end

% save path to repo
S.misc.main_path = pathRepo;

% do not run as batch job (parallel computing toolbox)
S.solver.run_as_batch_job = false;

% parallel cluster
S.solver.par_cluster_name = [];

% batch job additional paths
S.solver.batch_job_paths = {};

%%
if ~isempty(varargin)

    reference_path = fullfile(pathRepo,'Subjects',varargin{1},['settings_',varargin{1},'.m']);

    if isfile(reference_path)
        disp(['Initialising settings from "',reference_path,'".'])
        run(reference_path);
    else
        warning(['Could not initialise from "',reference_path,'". ',...
            'Ignoring input argument "',varargin{1},'".']);
    end


end


end
