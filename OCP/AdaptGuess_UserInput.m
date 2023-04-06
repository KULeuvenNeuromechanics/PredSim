function [guess,bounds] = AdaptGuess_UserInput(S,guess,bounds)
% --------------------------------------------------------------------------
% AdaptGuess_UserInput
%   This function adapts the initial guess of muscle activations to comply
%   with the bounds.
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - guess -
%   * initial guess values for all optimisation variables
% 
%   - bounds -
%   * boundaries for all optimisation variables
%
% OUTPUT:
%   - guess -
%   * initial guess values for all optimisation variables
% 
%   - bounds -
%   * boundaries for all optimisation variables
% 
% Original author: Maarten Afschrift
% Original date: 21/Oct/2020
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

for i=1:length(bounds.a.lower)
    guess.a(guess.a(:,i) < bounds.a.lower(i),i) = bounds.a.lower(i);
end

% if isempty(bounds.tf.lower)
%     warning('No lower bound on final time defined. Using initial guess as lower bound')
%     bounds.tf.lower = guess.tf;
% end
% if isempty(bounds.tf.upper)
%     warning('No upper bound on final time defined. Using initial guess as upper bound')
%     bounds.tf.upper = guess.tf;
% end

end
