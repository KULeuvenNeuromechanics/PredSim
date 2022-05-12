function [massM,tensions] = GetMuscleMass(FMo,lMo,tensions)
% --------------------------------------------------------------------------
% GetMuscleMass
%   compute the mass of the muscle
%   
% INPUT:
%   - FMo -
%   * 
%
%   - lMo -
%   * 
% 
%   - tensions -
%   * 
%
% OUTPUT:
%   - massM -
%   * 
%
%   - tensions -
%   * 
% 
% Original author: 
% Original date: 
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

volM = FMo.*lMo;
massM = volM.*(1059.7)./(tensions*1e6);

end

