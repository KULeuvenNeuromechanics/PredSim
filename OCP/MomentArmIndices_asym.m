% This function returns the indices of the muscles actuating the different 
% joints (for use with the moment arms) and the indices of the muscles in 
% the vector containing all muscle names.
%
% Author: Antoine Falisse
% Modified by Friedl De Groote for asymmetric model
% Date: 12/19/2018
% 
function [Indmusi,mai] = ...
    MomentArmIndices_asym(muscleNames,muscle_spanning_joint_INFO)

for i = 1:length(muscleNames)
    % Muscles are ordered as left first and right second
    Indmusi.(muscleNames{i}(1:end-2)) = find(strcmp(muscleNames,muscleNames{i}));
    % We select the muscles that are actuating the joints
    % In muscle_spanning_joint_INFO: the columns are organized as follows:
    % 1) hip flex r, 2) hip add r, 3) hip rot r, 4) knee flex r, 5) ankle flex r, 
    % 6) sub r, 7) hip flex l, 8) hip add l, 9) hip rot l, 10) knee flex l, 
    % 11) ankle flex l, 12) sub l, 13) trunk ext, 14) trunk ben, 15) trunk rot
    for j = 1:size(muscle_spanning_joint_INFO,2)
        if (muscle_spanning_joint_INFO(i,j) == 1)
           mai(j).mus(1,i) = Indmusi.(muscleNames{i}(1:end-2));
        end        
    end    
end

for j = 1:size(muscle_spanning_joint_INFO,2)
    mai(j).mus(mai(j).mus == 0) = [];
end

end
