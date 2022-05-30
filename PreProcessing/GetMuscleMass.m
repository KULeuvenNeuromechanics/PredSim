function [massM,specific_tension] = GetMuscleMass(FMo,lMo,specific_tension)
% --------------------------------------------------------------------------
% GetMuscleMass
%   This function computes the mass of a muscle.
%   
% INPUT:
%   - FMo -
%   * Optimal active muscle force
%
%   - lMo -
%   * Optimal muscle fiber length
% 
%   - specific_tension -
%   * Specific tension of muscle fibers
%
% OUTPUT:
%   - massM -
%   * mass of a muscle
%
%   - specific_tension -
%   * Specific tension of muscle fibers
% 
% Original author: Dhruv Gupta
% Original date: 
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

volM = FMo.*lMo;
massM = volM.*(1059.7)./(specific_tension*1e6);

end

