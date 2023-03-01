function hash = get_git_hash
% --------------------------------------------------------------------------
%get_git_hash 
%   gets the hash of the current github code used
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
% Last edit by: /
% Last edit date: DD/MM/YYYY
% --------------------------------------------------------------------------

command = [ 'git rev-parse HEAD' ];
[status,hash] = system(command);

if( status ~= 0 )
    error('Unable to get hash from file.');
end

end