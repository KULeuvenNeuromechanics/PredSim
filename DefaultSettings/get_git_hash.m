function [local_hash, branch_name, remote_hash] = get_git_hash(pathRepo)
% --------------------------------------------------------------------------
%get_git_hash 
%   Gets the hash of the current github code used. The git hash is used to
%   identify a specific version of the code. You can find the specific
%   version's code on GitHub by pasting the hash at 'commit-hash':
%   https://github.com/KULeuvenNeuromechanics/PredSim/tree/commit-hash
% 
% INPUT:
%   -pathRepo-
%   * path to the repository
% 
% OUTPUT:
%   -local_hash-
%   * the full hash for the commit you are currently using
%
%   -branch_name-
%   * the branch you are currently using
%
%   -remote_hash-
%   * the full hash for the latest commit on the remote
% 
% Original author: Bram Van Den Bosch
% Original date: 28/02/2023
%
% Last edit by: Bram Van Den Bosch
% Last edit date: 26/05/2023
% --------------------------------------------------------------------------

% save current working directory
pathInitDir = pwd;
% move to PredSim
cd(pathRepo)

% get hash of the local instance
command = ['git rev-parse HEAD'];
[status,local_hash] = system(command);

if contains(local_hash,"'git' is not recognized as an internal or external command")
    warning('Unable to get git hash. Git seems not to be installed on your machine or cannot be executed from the command line.');
elseif status ~= 0 
    warning('Unable to get git hash. It is advised to get PredSim through GitHub to have version control and to receive future updates.');
end

if status == 0
    % get name of the current branch
    command = ['git branch --show-current'];
    [~,branch_name] = system(command);
    
    % get the hash of the latest commit on the remote
    command = ['git rev-parse origin/' branch_name];
    [~,remote_hash] = system(command);
    
    % display warning, pointing to the most recent commit on the remote
    if ~strcmp(local_hash,remote_hash)
        warning(['There is a more recent version of this branch on the remote: https://github.com/KULeuvenNeuromechanics/PredSim/tree/' remote_hash]);
    end
end

% return to initial directory
cd(pathInitDir)

end