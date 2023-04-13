function hash = get_git_hash
% --------------------------------------------------------------------------
%get_git_hash 
%   Gets the hash of the current github code used. The git hash is used to
%   identify a specific version of the code. You can find the specific
%   version's code on GitHub by pasting the hash at 'commit-hash':
%   https://github.com/KULeuvenNeuromechanics/PredSim/tree/commit-hash
% 
% INPUT: none
% 
% OUTPUT:
%   -hash-
%   * the full hash for the commit you are currently using
% 
% Original author: Bram Van Den Bosch
% Original date: 28/02/2023
%
% Last edit by: Bram Van Den Bosch
% Last edit date: 13/04/2023
% --------------------------------------------------------------------------

command = ['git rev-parse HEAD'];
[status,hash] = system(command);

if contains(hash,"'git' is not recognized as an internal or external command")
    warning('Unable to get git hash. Git seems not to be installed on your machine or cannot be executed from the command line.')
elseif status ~= 0 
    warning('Unable to get git hash. It is advised to get PredSim through GitHub to have version control and to receive future updates.');
end

end