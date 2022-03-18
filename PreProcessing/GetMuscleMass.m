function [massM,tensions] = GetMuscleMass(FMo,lMo,tensions)

% compute the mass of the muscle
volM = FMo.*lMo;
massM = volM.*(1059.7)./(tensions'*1e6);

end

