clear all
close all

%%

osim_file = "C:\GBW_MyPrograms\PredSim_2022\AdaptOpenSimModel\subject1_mtpPin_unlocked_os4.osim";

%% initialize model

import org.opensim.modeling.*

%Create the Original OpenSim model from a .osim file

filename = osim_file;
ModelIn = Model(filename);
% StateIn = ModelIn.initSystem;

%%
osimFilexmlDoc = xmlread(osim_file);

nM = osimFilexmlDoc.getElementsByTagName('ForceSet').item(0).getElementsByTagName('objects').item(0).getElementsByTagName('Thelen2003Muscle').getLength;

% maybe first find indices of the quad muscles we want to adept? And then
% use the loop below to loop through them.

sp_muscles = {'rect_fem_r','vas_med_r','vas_int_r','vas_lat_r',...
    'rect_fem_l','vas_med_l','vas_int_l','vas_lat_l'}; % muscles with simmspline
sp_parents = {'tibia_r','tibia_r','tibia_r','tibia_r',...
    'tibia_l','tibia_l','tibia_l','tibia_l'};
sp_joints = {'knee_r','knee_r','knee_r','knee_r',...
    'knee_l','knee_l','knee_l','knee_l'};

order = 4;
a = -2.5:0.01:2.5;

for i=1:nM
    for j = 1:length(sp_muscles)
        if strcmp(char(osimFilexmlDoc.getElementsByTagName('Thelen2003Muscle').item(i-1).getAttributes.item(0).getValue),sp_muscles(j))
            mp_x_loc = osimFilexmlDoc.getElementsByTagName('Thelen2003Muscle').item(i-1).getElementsByTagName('x_location').item(0); % muscle point x location
            mp_y_loc = osimFilexmlDoc.getElementsByTagName('Thelen2003Muscle').item(i-1).getElementsByTagName('y_location').item(0); % muscle point y location
            mp_z_loc = osimFilexmlDoc.getElementsByTagName('Thelen2003Muscle').item(i-1).getElementsByTagName('z_location').item(0); % muscle point z location
    
            spline(j).x_location.x = str2num(mp_x_loc.getElementsByTagName('SimmSpline').item(0).getElementsByTagName('x').item(0).getFirstChild.getNodeValue);
            spline(j).x_location.y = str2num(mp_x_loc.getElementsByTagName('SimmSpline').item(0).getElementsByTagName('y').item(0).getFirstChild.getNodeValue);
            spline(j).x_location.s = str2num(mp_x_loc.getElementsByTagName('scale').item(0).getFirstChild.getNodeValue);
            
            spline(j).y_location.x = str2num(mp_y_loc.getElementsByTagName('SimmSpline').item(0).getElementsByTagName('x').item(0).getFirstChild.getNodeValue);
            spline(j).y_location.y = str2num(mp_y_loc.getElementsByTagName('SimmSpline').item(0).getElementsByTagName('y').item(0).getFirstChild.getNodeValue);         
            spline(j).y_location.s = str2num(mp_y_loc.getElementsByTagName('scale').item(0).getFirstChild.getNodeValue);

            spline(j).z_location.x = str2num(mp_z_loc.getElementsByTagName('SimmSpline').item(0).getElementsByTagName('x').item(0).getFirstChild.getNodeValue);
            spline(j).z_location.y = str2num(mp_z_loc.getElementsByTagName('SimmSpline').item(0).getElementsByTagName('y').item(0).getFirstChild.getNodeValue);           
            spline(j).z_location.s = str2num(mp_z_loc.getElementsByTagName('scale').item(0).getFirstChild.getNodeValue);

            
            spline(j).poly.x = polyfit(spline(j).x_location.x,spline(j).x_location.y*spline(j).x_location.s,order);
            spline(j).poly.y = polyfit(spline(j).y_location.x,spline(j).y_location.y*spline(j).y_location.s,order);
            spline(j).poly.z = polyfit(spline(j).z_location.x,spline(j).z_location.y*spline(j).z_location.s,0);

        end
    end
end

%%
for i = 1:length(sp_muscles)
    % find pathpoints (future work: automatically search for pathpoints
    % with simmspline)
    name_currentmuscle = sp_muscles{1,i};
    currentMuscle = ModelIn.getMuscles.get(name_currentmuscle);
    pathpoints = currentMuscle.getGeometryPath.getPathPointSet;
    n_pp = pathpoints.getSize; % number of path points in osim model

    % remove pathpoint with simmspline (in the model I use now the last
    % pathpoint is always the one with the SimmSplines)
    pp_a =  pathpoints.get(n_pp-1);
    pathpoints.remove(pp_a);
    ModelIn.updForceSet;

    % add movingpathpoint with polynomial
    name   = [name_currentmuscle '-P' num2str(n_pp)];
    parent = sp_parents(i);
    joint  = sp_joints(i);
    poly   = spline(i).poly;
    addMovingPathPoint(ModelIn,name_currentmuscle,name,parent,joint,poly);
end
%%
ModelIn.finalizeConnections;
ModelIn.print("C:\GBW_MyPrograms\PredSim_2022\AdaptOpenSimModel\subject1_mtpPin_unlocked_os4_test.osim");
