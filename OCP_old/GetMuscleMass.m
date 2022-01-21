function [massM,tensions] = GetMuscleMass(muscleNames,params)
% get the muscle parameters
FMo = params(1,:);
lMo = params(2,:);

% get the specific tension for each muscle
tensions = getSpecificTensions(muscleNames);

% compute the mass of the muscle
volM = FMo.*lMo;
massM = volM.*(1059.7)./(tensions*1e6);

end

