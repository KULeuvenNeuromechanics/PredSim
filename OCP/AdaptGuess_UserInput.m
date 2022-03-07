function [guess,bounds] = AdaptGuess_UserInput(guess,bounds,S)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
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
