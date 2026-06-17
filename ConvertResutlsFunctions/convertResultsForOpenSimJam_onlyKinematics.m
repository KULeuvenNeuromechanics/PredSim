% take the outputs from pred sim and pull them out for use in opensim
% specifically I want to pull out muscle activations, forces, and ground
% reactions for re-running kineatics, and verification of joint contact
% loading reductions


% typcial
clear all
close all
clc
%% what do I need to do .

folder_dll = 'C:\Users\u0130218\Documents\MATLAB\OpenSimJAM\opensim-core-4.3-2021-06-27-54b40380c\bin';

OLD_path = getenv('PATH');
if isempty(strfind(getenv('PATH'), folder_dll))
    setenv('PATH', [folder_dll ';' getenv('PATH')]);
end

import org.opensim.modeling.*

% define the base dir
kinMotDir = 'C:\Users\u0130218\Documents\MATLAB\predsim_runInOpenSimJam\ISB\AverageKinematics_S204_R.mot';

% start by laoding the kinemaitcs file to get teh time.
[kinMotData] = read_motionFile_v40(kinMotDir);

lookupTable = struct();
lookupTable.('predsimToJAM') = struct( 'pelvis_rotation','pelvis_rot' , 'hip_flexion_r','hip_flex_r',...
    'hip_flexion_l','hip_flex_l' , 'hip_adduction_r','hip_add_r' , 'hip_adduction_l','hip_add_l',...
    'hip_rotation_r','hip_rot_r' , 'hip_rotation_l','hip_rot_l',...
    'knee_angle_r','knee_flex_r' , 'knee_angle_l','knee_flex_l' , ...
    'ankle_angle_r','ankle_flex_r', 'ankle_angle_l', 'ankle_flex_l',...
     'lumbar_extension', 'lumbar_ext' ,'lumbar_bending','lumbar_latbend','lumbar_rotation','lumbar_rot',...
    'subtalar_angle_l', 'subt_angle_l','subtalar_angle_r', 'subt_angle_r');

% return the list of predsim labels to change
labelsToUpd = fieldnames(lookupTable.predsimToJAM);
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
        kinMotData.labels(fl) = {lookupTable.predsimToJAM.(char(oldLabel))}; 
        end     
    end      
end

toRev = {'knee_flex_r' , 'knee_flex_l'};

for x = 1:length(toRev)
    tmptoRev = toRev(x);
    for l = 1:length(kinMotData.labels)
       tmpheader= kinMotData.labels(l);
       if strcmp(tmptoRev,tmpheader)
           kinMotData.data(:,l)= kinMotData.data(:,l)*-1;
       end   
    end   
end
% add the initial guess for second coords
% load the spline functions

% uses a dummy file right now but can be personalised 
splineFunctions = FunctionSet('secondary_coordinate_constraint_functions_r.xml');

for ci = 0:splineFunctions.getSize()-1
    tempFunc = SimmSpline.safeDownCast(splineFunctions.get(ci));
    tempFuncName = char(tempFunc.getName());
    tempFuncCoordName = strsplit(tempFuncName,'/');
    tempFuncCoordName = char(tempFuncCoordName(4));
    % add it to the kin mot dir
    kinMotData.labels(end+1) = {tempFuncCoordName};

    % loop through the knee fe 
    
    knee_index = find(contains(fullLabels,'knee_angle_r'));
    kneeVals = kinMotData.data(:,knee_index);
    tvals = nan(length(kneeVals),1);
    for kci  = 1:length(kneeVals)
    %tempFunc.calcValue(Vector(1,deg2rad(tempDeg))))
     tvals(kci,1) = tempFunc.calcValue(Vector(1,deg2rad(kneeVals(kci))));
    end
    
    % append the data
    kinMotData.data(:,end+1) = tvals;
end

% uses a dummy file right now but can be personalised 
splineFunctions = FunctionSet('secondary_coordinate_constraint_functions_l.xml');

for ci = 0:splineFunctions.getSize()-1
    tempFunc = SimmSpline.safeDownCast(splineFunctions.get(ci));
    tempFuncName = char(tempFunc.getName());
    tempFuncCoordName = strsplit(tempFuncName,'/');
    tempFuncCoordName = char(tempFuncCoordName(4));
    % add it to the kin mot dir
    kinMotData.labels(end+1) = {tempFuncCoordName};

    % loop through the knee fe 
    
    knee_index = find(contains(fullLabels,'knee_angle_l'));
    kneeVals = kinMotData.data(:,knee_index);
    tvals = nan(length(kneeVals),1);
    for kci  = 1:length(kneeVals)
    %tempFunc.calcValue(Vector(1,deg2rad(tempDeg))))
     tvals(kci,1) = tempFunc.calcValue(Vector(1,deg2rad(kneeVals(kci))));
    end
    
    % append the data
    kinMotData.data(:,end+1) = tvals;
end
% write out
write_motionFile_v40(kinMotData, strrep(kinMotDir, '.mot', '_JAM.mot'));

    