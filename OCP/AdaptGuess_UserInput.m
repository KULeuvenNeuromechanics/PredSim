function [guess,bounds] = AdaptGuess_UserInput(S,guess,bounds)
% --------------------------------------------------------------------------
% AdaptGuess_UserInput
%   (Explanation)
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - guess -
%   * 
% 
%   - bounds -
%   * 
%
% OUTPUT:
%   - guess -
%   * 
% 
%   - bounds -
%   * 
% 
% Original author: 
% Original date: 
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

guess.a(guess.a < S.bounds.a.lower) = S.bounds.a.lower;
if isempty(bounds.tf.lower)
    warning('No lower bound on final time defined. Using initial guess as lower bound')
    bounds.tf.lower = guess.tf;
end
if isempty(bounds.tf.upper)
    warning('No upper bound on final time defined. Using initial guess as upper bound')
    bounds.tf.upper = guess.tf;
end
end
