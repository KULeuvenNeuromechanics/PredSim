% This function returns the indices of the muscles actuating the different 
% joints (for use with the moment arms) and the indices of the muscles in 
% the vector containing all muscle names.
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function [Indmusi,mai] = ...
    MomentArmIndices(muscleNames,muscle_spanning_joint_INFO)

DOFnames = {'hip flex','Hip add','hip rot','knee flex','ankle flex','subtalar',...
    'mtp angle','trunk extension','trunk bending','trunk rotation'};

NMuscle = length(muscleNames)*2;
for i = 1:length(muscleNames)
    % Muscles are ordered as left first and right second
    Indmusi.(muscleNames{i}(1:end-2)).l = find(strcmp(muscleNames,muscleNames{i}));
    Indmusi.(muscleNames{i}(1:end-2)).r = Indmusi.(muscleNames{i}(1:end-2)).l + NMuscle/2;
    % We select the muscles that are actuating the joints
    % In muscle_spanning_joint_INFO: the columns are organized as follows:
    % 1) hip flex, 2) hip add, 3) hip rot, 4) knee flex, 5) ankle flex, 
    % 6) sub, 7) trunk ext, 8) trunk ben, 9) trunk rot
    for j = 1:size(muscle_spanning_joint_INFO,2)
        if (muscle_spanning_joint_INFO(i,j) == 1)
           mai(j).mus.l(1,i) = Indmusi.(muscleNames{i}(1:end-2)).l;
           mai(j).mus.r(1,i) = Indmusi.(muscleNames{i}(1:end-2)).r;
        end
        if i ==1
            mai(j).dof = DOFnames{j};
        end
    end    
end

for j = 1:size(muscle_spanning_joint_INFO,2)
    mai(j).mus.l(mai(j).mus.l == 0) = [];
    mai(j).mus.r(mai(j).mus.r == 0) = [];
end

end
