function [massM,tensions] = GetMuscleMass(muscleNames,params,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
BoolRajagopal = 0;
if ~isempty(varargin)
    BoolRajagopal = varargin{1};
end

% get the muscle parameters
FMo = params(1,:);
lMo = params(2,:);

% get the specific tension for each muscle
if BoolRajagopal
    tension = getSpecificTensions(muscleNames);
else
    tension = getSpecificTensions(muscleNames(1:end-3));
end
tensions = [tension;tension]';

% compute the mass of the muscle
volM = FMo.*lMo;
massM = volM.*(1059.7)./(tensions*1e6);

end

