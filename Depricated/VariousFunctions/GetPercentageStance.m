function [PartStance, PartSwing, DoubleSupport ] = GetPercentageStance(Fy,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if~isempty(varargin)
    threshold = varargin{1};
else
    threshold = 30;
end

% get stance and swing phase
nFR = length(Fy(:,1));
PartStance =  sum(Fy(:,1)>threshold)./nFR*100;
PartSwing = 100-PartStance;

% get the double support phase
DoubleSupport = sum(Fy(:,1)>threshold & Fy(:,2)>threshold)./nFR*100;



end

