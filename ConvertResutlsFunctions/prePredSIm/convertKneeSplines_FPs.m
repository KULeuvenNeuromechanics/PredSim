

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


Model = Model('C:\Users\u0130218\Documents\Data\FPmodels\combinedModel.osim');
%se = dummyJamModel.getSimbodyEngine();
State = Model.initSystem();

% load the motion file 
%
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop through frames
%                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kinData = read_motionFile_v40('C:\Users\u0130218\Documents\MATLAB\predsim_runInOpenSimJam\postESMAC\woTracking\med10lat10\jamOutput\comak\walking_kinematics.sto');
% strucutre = kinData.labels ; kinData.data


% create dummy file for testing
kinData= struct();
kinData.labels = {'time', '/jointset/osimjam_knee_r/osimjam_knee_flex_r/value','/jointset/osimjam_knee_r/osimjam_knee_add_r/value','/jointset/osimjam_knee_r/osimjam_knee_rot_r/value' ... 
    '/jointset/osimjam_knee_r/osimjam_knee_tx_r/value','/jointset/osimjam_knee_r/osimjam_knee_ty_r/value','/jointset/osimjam_knee_r/osimjam_knee_tz_r/value' ...
    '/jointset/predsim_knee_r/predsim_knee_angle_r/value', '/jointset/predsim_knee_r/predsim_knee_add_r/value', '/jointset/predsim_knee_r/predsim_knee_rot_r/value'... 
    '/jointset/predsim_knee_r/predsim_knee_tx_r/value', '/jointset/predsim_knee_r/predsim_knee_ty_r/value', '/jointset/predsim_knee_r/predsim_knee_tz_r/value'
    };

%kinData.data = [0, -90,10,25,...
%                90,10,25, ...
%                0.01, 0.02, 0.005];

kneeflexions =deg2rad(-10:1:120);
nr = length(kneeflexions);
times = 0:0.05:(0.05*(nr-1));
kinData.inDeg='no';

% add times to data
kinData.data(:,1) = times';
% add flexions for osimjam
kinData.data(:,2) = kneeflexions';
% add flexions for predsim
kinData.data(:,8) = kneeflexions*-1';

% create a lookup table for the DOFS
lut=struct();
lut.('assignCoords') = struct('knee_add_r',[3,9],'knee_rot_r',[4,10],'knee_tx_r',[5],'knee_ty_r',[6],'knee_tz_r',[7]);


% uses a dummy file right now but can be personalised 
splineFunctions = FunctionSet('secondary_coordinate_constraint_functions_r.xml');

 
%for ti = 1:kinData.nr
%kinData.inDeg = 'no';
    
for ci = 0:splineFunctions.getSize()-1
    tempFunc = SimmSpline.safeDownCast(splineFunctions.get(ci));
    tempFuncName = char(tempFunc.getName());
    tempFuncCoordName = strsplit(tempFuncName,'/');
    tempFuncCoordName = char(tempFuncCoordName(4));
    % add it to the kin mot dir

    for kci  = 1:length(kneeflexions)
        tvals(kci,1) = tempFunc.calcValue(Vector(1,kneeflexions(kci)));
    end
    % secondary coord index
    try
        tsci = lut.('assignCoords').(tempFuncCoordName);
        % append the data
        %kinMotData.data(:,end+1) = tvals;
        for i = 1:length(tsci)
            kinData.data(:,tsci(i)) = tvals;
        end
    catch
       disp(append(' Skipping the following dof ..' , char(tempFuncCoordName) ))
    end
end
%kinData.labels = kinData.labels(1:10);
%write_motionFile_v40(kinData, 'test.mot');
    
%%% NOW LOOP THROUGH ADN SET ALL THE LABELS 
kinToSet = kinData.labels(2:10);   

for kfr = 1:length(kinData.data);
    
    for kts = 1:length(kinToSet)
        % cut the crap off the end to get the coord
        coordName = strsplit(char(kinToSet(kts)),'/');
        coordName = char(coordName(4));
        % get the index in data
        Tindex = find(contains(kinData.labels,kinToSet(kts)));
        Tval = kinData.data(kfr,Tindex);
        % set the model pose
        Model.getCoordinateSet().get(coordName).setValue(State, Tval);
        State = Model.updWorkingState();
    end
    
        % extract the frames 
        %parFrame = Model.getBodySet().get('femur_distal_r').getTransformInGround(State);
        predsim_chiFrame = Model.getBodySet().get('predsim_tibia_proximal_r').getTransformInGround(State);
        osimjam_chiFrame = Model.getBodySet().get('osimjam_tibia_proximal_r').getTransformInGround(State);
        
        % calculate the roations 
 
       % parRotO = parFrame.R();
       % parTransO = parFrame.T();
        % conv to mlab array
        %parRot = convOsimMatToMat(parRotO); 
        
        % get the values for the child - tibia proximal
        
        predsim_chiRotO = predsim_chiFrame.R();
        predsim_chiTransO = predsim_chiFrame.T();
        % conv to mlab array
        %predsim_chiRot = convOsimMatToMat(predsim_chiRotO);   

        osimjam_chiRotO = osimjam_chiFrame.R();
        osimjam_chiTransO = osimjam_chiFrame.T();
        % conv to mlab array
        %osimjam_chiRot = convOsimMatToMat(osimjam_chiRotO);   
        
        for xx = 1:3
            predsim_chiTrans(xx)=predsim_chiTransO.get(xx-1);
            osimjam_chiTrans(xx)= osimjam_chiTransO.get(xx-1);
        end
        
        % calcualte the difference
        transCorr = osimjam_chiTrans - predsim_chiTrans;
        % 
        kinData.data(kfr,11) = transCorr(1);
        kinData.data(kfr,12) = transCorr(2);
        kinData.data(kfr,13) = transCorr(3);
end
%end frame loop

write_motionFile_v40(kinData, 'ConvertedData.mot');


%% Now this works... 
% convert the results from kinemaitcs to polynomials 

% in theory I can hard code this 
indepInd= 8;
% now deps
depInds = [9,10,11,12,13];

% graph poly fit values
gr=0;

for di = 1:length(depInds)
   
    % pull the data 
    indData = kinData.data(:,indepInd);
    depData = kinData.data(:,depInds(di));
   
    tempFuncCoordName = strsplit(char(kinData.labels(depInds(di))),'/');
    tempFuncCoordName = char(tempFuncCoordName(4));
    
    % poly fit this bitch
    p = polyfit(indData , depData, 4);
    if gr == 1
        h = figure();
        ny = polyval(p, indData);
        plot(indData, depData,'r')
        hold on
        plot(indData, ny,'b')
        waitfor(h)
    end
    
    fprintf(' Polyvals for %s are %1.5f %1.5f %1.5f %1.5f %1.5f, \n', tempFuncCoordName,  p(1),p(2), p(3), p(4), p(5))
end





