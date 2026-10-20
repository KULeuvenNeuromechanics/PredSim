% take the outputs from pred sim and pull them out for use in opensim
% specifically I want to pull out muscle activations, forces, and ground
% reactions for re-running kineatics, and verification of joint contact
% loading reductions


% typcial
clear all
close all
clc
%% what do I need to do .

kinMotDir = 'C:/Users/u0130218/Documents/Data/HappyJoints_PredSimTesting/S204/COMAK/Gait/Right/D/comak-inverse-kinematics/ik_generic/S204_ik.mot';
% if it is a right leg trial - you need to conv from Rad the left knee 
toConv = {'knee_angle_l'};

kinOutMotDir = strrep(kinMotDir, '.mot', '_forPredSim.mot');
% start by laoding the kinemaitcs file to get teh time.
[kinMotData] = read_motionFile_v40(kinMotDir);

% create a strcut to convert.

lookupTable = struct();
lookupTable.('JAMtoPredSim') = struct( 'pelvis_rot' ,'pelvis_rotation', 'hip_flex_r','hip_flexion_r',...
    'hip_flex_l' ,'hip_flexion_l', 'hip_add_r' ,'hip_adduction_r', 'hip_add_l','hip_adduction_l',...
    'hip_rot_r' ,'hip_rotation_r', 'hip_rot_l','hip_rotation_l',...
    'knee_flex_r' ,'knee_angle_r', 'knee_flex_l' ,'knee_angle_l', ...
    'ankle_flex_r','ankle_angle_r', 'ankle_flex_l', 'ankle_angle_l',...
    'lumbar_ext', 'lumbar_extension' ,'lumbar_latbend','lumbar_bending','lumbar_rot','lumbar_rotation',...
    'subt_angle_l', 'subtalar_angle_l','subt_angle_r', 'subtalar_angle_r');



% return the list of labels to change
labelsToUpd = fieldnames(lookupTable.JAMtoPredSim);
% full
fullLabels = kinMotData.labels;
% loop throug hthisese
for ol = 1:length(labelsToUpd)
    % loop through and check
    oldLabel = labelsToUpd(ol);
    for fl =1:length(fullLabels)
        tmpLabel = fullLabels(fl);
        if strcmp(oldLabel,tmpLabel)
        % if match upd
        kinMotData.labels(fl) = {lookupTable.JAMtoPredSim.(char(oldLabel))}; 
        end
        
    end   
    
end

toRev = {'knee_angle_r' , 'knee_angle_l', 'hip_adduction_l' , 'hip_rotation_l'};

for x = 1:length(toRev)
    tmptoRev = toRev(x);
    for l = 1:length(kinMotData.labels)
       tmpheader= kinMotData.labels(l);
       if strcmp(tmptoRev,tmpheader)
           kinMotData.data(:,l)= kinMotData.data(:,l)*-1;
       end    
    end
    
end


for x = 1:length(toConv)
    tmptoRev = toConv(x);
    for l = 1:length(kinMotData.labels)
       tmpheader= kinMotData.labels(l);
       if strcmp(tmptoRev,tmpheader)
           kinMotData.data(:,l)= rad2deg(kinMotData.data(:,l));
       end    
    end
    
end


% write out
write_motionFile_v40(kinMotData, kinOutMotDir);

disp('done cuh')
