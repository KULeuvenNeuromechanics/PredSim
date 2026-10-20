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
%
% --------------------------------------------------------------------------
% This file is part of PredSim.
% 
% PredSim: A Framework for Rapid Predictive Simulations of Locomotion
% Copyright (c) 2026 KU Leuven
% 
% PredSim is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Affero General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% 
% PredSim is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public 
% License for more details.
% 
% You should have received a copy of the GNU Affero General Public License 
% along with PredSim. If not, see <https://www.gnu.org/licenses/>.
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
