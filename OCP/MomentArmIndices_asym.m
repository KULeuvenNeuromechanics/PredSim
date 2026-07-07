function [Indmusi,mai] = MomentArmIndices_asym(muscleNames,muscle_spanning_joint_INFO)
% --------------------------------------------------------------------------
% MomentArmIndices_asym
%   This function returns the indices of the muscles actuating the different 
%   joints (for use with the moment arms) and the indices of the muscles in 
%   the vector containing all muscle names.
%   
% INPUT:
%   - muscleNames -
%   * cell array with names of muscles
%
%   - muscle_spanning_joint_INFO -
%   * table with a column for each coordinate and a row for each muscle. 1
%   means this muscle and coordinate interact, 0 means they don't
%
% OUTPUT:
%   - Indmusi -
%   * struct with muscle indices that cross each joint
%
%   - mai -
%   * indices of the coordinates that each muscle crosses
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
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

for j = 1:size(muscle_spanning_joint_INFO,2)
    c=0;
    for i = 1:length(muscleNames)
        Indmusi.(muscleNames{i}) = find(strcmp(muscleNames,muscleNames{i}));
        if (muscle_spanning_joint_INFO(i,j) == 1)
            c=c+1;
            mai(j).mus(1,c) = Indmusi.(muscleNames{i});
        end        
    end    
end

end
