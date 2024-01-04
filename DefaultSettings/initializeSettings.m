function [S] = initializeSettings(varargin)
% --------------------------------------------------------------------------
%initializeSettings
%   This function creates the empty settings struct S up to the field 
%   above the field containing data. 
% 
% INPUT:
%   - reference_name - (optional input)
%   * pass the name of a reference model to initialise the settings
%   according to the relevant publication. 
% 
% OUTPUT:
%   - S -
%   * empty settings struct S
%
% Original author: Bram Van Den Bosch
% Original date: 01/12/2021
%
% Last edit by: Bram Van Den Bosch
% Last edit date: 26/May/2023
% --------------------------------------------------------------------------

S = struct;

S.metabolicE   = [];
S.misc         = [];
S.post_process = [];
S.solver       = [];
S.subject      = [];
S.weights      = [];
S.Cpp2Dll      = [];
S.orthosis     = [];

% bounds have an .upper and .lower field
S.bounds.a          = [];
S.bounds.calcn_dist = [];
S.bounds.femur_hand_dist = [];
S.bounds.toes_dist  = [];
S.bounds.tibia_dist = [];
S.bounds.SLL        = [];
S.bounds.SLR        = [];
S.bounds.dist_trav  = [];
S.bounds.t_final    = [];

% polynomial order has .lower and .upper field
S.misc.poly_order = [];

% save the git hash
[S.misc.git.local_hash,S.misc.git.branch_name, S.misc.git.remote_hash] = get_git_hash;

% initiate for warning
S.subject.adapt_IG_pelvis_y = 0;


%%
if ~isempty(varargin)

    [pathDefaultSettings,~,~] = fileparts(mfilename('fullpath'));
    [pathRepo,~,~] = fileparts(pathDefaultSettings);

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