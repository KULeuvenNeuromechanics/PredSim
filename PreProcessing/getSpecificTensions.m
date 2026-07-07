function specific_tension = getSpecificTensions(S, muscleNames)
% --------------------------------------------------------------------------
% getSpecificTensions
%   Returns the specific tension value for each muscle.
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - muscleNames -
%   * Cell array of muscle names
%
% OUTPUT:
%   - specific_tension -
%   * Array with specific tension of each muscle. Default is 0.70
% 
% Original author: Lars D'Hondt and Tim van der Zee
% Original date: 16 September 2025
%
% --------------------------------------------------------------------------
% This file is part of PredSim.
% 
% PredSim: A Framework for Rapid Predictive Simulations of Locomotion
% Copyright (c) 2026 KU Leuven
% 
% PredSim is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Affero General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% 
% PredSim is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public 
% License for more details.
% 
% You should have received a copy of the GNU Affero General Public License 
% along with PredSim. If not, see <https://www.gnu.org/licenses/>.
% --------------------------------------------------------------------------

specific_tension = 0.25*ones(length(muscleNames),1);

if isnumeric(S.subject.default_specific_tension) && S.subject.default_specific_tension > 0
    specific_tension(:) = S.subject.default_specific_tension;

elseif exist(S.subject.default_specific_tension,'file')
    default_sigma = readtable(S.subject.default_specific_tension);

    for i=1:length(muscleNames)
        default_sigma_i = default_sigma(strcmp(default_sigma.name, muscleNames{i}),2);
        if ~isempty(default_sigma_i)
            specific_tension(i) = table2array(default_sigma_i);
        end
    end
end

end    
