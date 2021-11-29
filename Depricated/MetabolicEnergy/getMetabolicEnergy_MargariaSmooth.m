function [Metab] = getMetabolicEnergy_MargariaSmooth(Fce,vM,b)
%getMetabolicEnergy_MargariaSmooth Computes metabolic power of the muscles
%based on Margaria et al. 1968
%
%   Detailed explanation goes here
% [Margaria et al. 1968: Positive and negative work performances and their 
% efficiencies in human locomotion. Euro- pean journal of applied physiology
% and occupational physiology
% 
% Author: Maarten Afschrift


Power =  Fce .* -vM; % minus since shortening is negative

% get positive and negative power with tanh function (smooth approximation)
PowerPos = Power.*(tanh(b*Power)*0.5+0.5);      % negative power = 0, positive power = posivive power
PowerNeg = -Power.*(tanh(b*-Power)*0.5+0.5);    % negative power => positive value, postive power = 0;

% compute Metabolic power
Metab = PowerPos./0.25 + PowerNeg./1.2;







end

