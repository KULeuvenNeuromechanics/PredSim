function [Indmusi,mai] = ...
    MomentArmIndices_asym(muscleNames,muscle_spanning_joint_INFO)
% --------------------------------------------------------------------------
% MomentArmIndices_asym
%   This function returns the indices of the muscles actuating the different 
%   joints (for use with the moment arms) and the indices of the muscles in 
%   the vector containing all muscle names.
%   
% INPUT:
%   - muscleNames -
%   * 
%
%   - muscle_spanning_joint_INFO -
%   * 
%
% OUTPUT:
%   - Indmusi -
%   * 
%
%   - mai -
%   * 
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
%   Modified by Friedl De Groote for asymmetric model
%   Modified by Dhruv Gupta for efficiency
% Last edit by: 
% Last edit date: 
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
