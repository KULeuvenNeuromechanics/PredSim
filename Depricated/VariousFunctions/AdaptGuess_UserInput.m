function [guess] = AdaptGuess_UserInput(guess,bounds,S)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if ~isempty(S.Bounds)
    Inds = guess.a(1,:) < bounds.a.lower;
    for i=Inds
        guess.a(:,i) = bounds.a.lower(i);
    end
end
if ~isempty(S.Bounds.tf)
    guess.tf = nanmean(S.Bounds.tf);
end
end

