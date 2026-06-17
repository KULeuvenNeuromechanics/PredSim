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
baseResultsDir = 'C:\Users\u0130218\Documents\MATLAB\predsim_runInOpenSimJam\ISB\sims\effectTrackingPref\ineq\data';
% hard code 

itID = 'IG_PrevSol_S204_JCF2_Tracked_LR_HK5A20_OCPweight_25000_250';

kinMotDir = fullfile(baseResultsDir , strcat(itID,'.mot'));

outMatDir = fullfile(baseResultsDir , strcat(itID,'.mat'));

% load the results
outMat = load(outMatDir);

% just to verify
    
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
write_motionFile_v40(kinMotData, fullfile(baseResultsDir , strcat(itID,'_JAM.mot')));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start with the muslce activation/forces 

% muscle headers
muscHeaders = [{'time'},outMat.model_info.muscle_info.muscle_names];

% create the matrix values
actOutMatrix = [outMat.R.muscles.a;outMat.R.muscles.a];
% add time
actOutMatrix = [ kinMotData.data(:,1) , actOutMatrix];
% make it a struct for output
actOutData=struct();
actOutData.labels = muscHeaders; 
actOutData.data =  actOutMatrix; 
actOutData.inDeg = 'no'; 
% write out
write_motionFile_v40(actOutData, fullfile(baseResultsDir , strcat(itID,'_activations.sto')));


forceOutMatrix = [outMat.R.muscles.FT;outMat.R.muscles.FT];
% add time
forceOutMatrix = [ kinMotData.data(:,1) , forceOutMatrix];
% make it a struct for output
forceOutData=struct();
forceOutData.labels = muscHeaders; 
forceOutData.data =  forceOutMatrix; 
forceOutData.inDeg = 'no'; 
% write out
write_motionFile_v40(forceOutData, fullfile(baseResultsDir , strcat(itID,'_forces.sto')));

% now deal with the GRFs 
grfHeaders=[{'time'} , ...
    {'groundr_force_vx'},{'groundr_force_vy'},{'groundr_force_vz'}, ...
    {'groundr_force_px'},{'groundr_force_py'},{'groundr_force_pz'}, ...
    {'groundr_torque_x'},{'groundr_torque_y'},{'groundr_torque_z'}, ...
    {'groundl_force_vx'},{'groundl_force_vy'},{'groundl_force_vz'}, ...
    {'groundl_force_px'},{'groundl_force_py'},{'groundl_force_pz'},...
    {'groundl_torque_x'},{'groundl_torque_y'},{'groundl_torque_z'}
    ];

grfOutMatrix = [kinMotData.data(:,1), ...
    [outMat.R.ground_reaction.GRF_r;outMat.R.ground_reaction.GRF_r]...
    [outMat.R.ground_reaction.COP_r;(outMat.R.ground_reaction.COP_r + [kinMotData.data(101,5),0,0])]...
    [outMat.R.ground_reaction.GRM_r;outMat.R.ground_reaction.GRM_r]...
    [outMat.R.ground_reaction.GRF_l;outMat.R.ground_reaction.GRF_l]...
    [outMat.R.ground_reaction.COP_l;(outMat.R.ground_reaction.COP_l+ [kinMotData.data(101,5),0,0])]...
    [outMat.R.ground_reaction.GRM_l;outMat.R.ground_reaction.GRM_l]...
    ];

% adjsut GRF moments to be 0
grfOutMatrix(:,8) = 0; 
grfOutMatrix(:,9) = 0; 
grfOutMatrix(:,10) = 0; 

grfOutMatrix(:,17) = 0; 
grfOutMatrix(:,18) = 0; 
grfOutMatrix(:,19) = 0; 

% make it a struct for output
grfOutData = struct();
grfOutData.labels = grfHeaders;
grfOutData.data = grfOutMatrix;
grfOutData.inDeg = 'no'; 

% write out
write_motionFile_v40(grfOutData, fullfile(baseResultsDir , strcat(itID,'_grf.mot')));

% Output the inverse dyanmics  
idOutData = struct();
idOutData.labels = [{'time'}, outMat.R.colheaders.coordinates'];
idOutData.data = [kinMotData.data(:,1) ,[outMat.R.kinetics.T_ID;outMat.R.kinetics.T_ID]] ;
idOutData.inDeg = 'no';
% write out
write_motionFile_v40(idOutData, fullfile(baseResultsDir , strcat(itID,'_inverseDynamics.mot')));
    