function [] = playVideo(osim_file,mot_file,path_geom)

% kill visualizer before opening new one to prevent memory leaks
try
    [~,~] = system('taskkill /IM simbody-visualizer.exe');
catch
end

import org.opensim.modeling.*

osimModel = Model(osim_file);
osimModel.setUseVisualizer(true);
model_state = osimModel.initSystem();
vis = osimModel.getVisualizer();
vis.addDirToGeometrySearchPaths(path_geom);
pause(0.1);
vis.show(model_state);
% pause(3);

dat = read_motionFile_v40(mot_file);

dt = dat.data(2,1)-dat.data(1,1);
N = length(dat.data(:,1));                

% CoordNamesAPI = {'pelvis_tilt','pelvis_tx','pelvis_ty','hip_flexion_r','hip_flexion_l',...
%     'lumbar_extension','knee_angle_r','knee_angle_l','arm_flex_r','arm_flex_l',...
%     'ankle_angle_r','ankle_angle_l','elbow_flex_r','elbow_flex_l'};

if contains(osim_file, 'Vitruvian_Man.osim')
CoordNamesAPI = {'pelvis_tilt','pelvis_tx','pelvis_ty','hip_flexion_r','hip_flexion_l',...
    'lumbar_extension','knee_angle_r','knee_angle_l','arm_flex_r','arm_flex_l',...
    'ankle_angle_r','ankle_angle_l'};

else
CoordNamesAPI = {'pelvis_tilt','pelvis_tx','pelvis_ty','hip_flexion_r','hip_flexion_l',...
    'lumbar_extension','knee_angle_r','knee_angle_l',...
    'ankle_angle_r','ankle_angle_l'};
end

% column index in datastrcture for every coordName
IndexCoord = nan(length(CoordNamesAPI),1);
for i=1:length(CoordNamesAPI)
    IndexCoord(i) = find(strcmp(CoordNamesAPI{i},dat.labels))-1;
end
dat.data(:,IndexCoord(2)+1) = dat.data(:,IndexCoord(2)+1) - mean(dat.data(:,IndexCoord(2)+1));

simbodyVis = vis.updSimbodyVisualizer();
simbodyVis.setShowSimTime(true);

Qvect = model_state.getQ();
for i=[1:N 1]
    dSel = dat.data(i,2:end);
    for j=1:length(CoordNamesAPI)
        if j==2 || j==3
            Qvect.set(j-1,dSel(IndexCoord(j)));
        else
            Qvect.set(j-1,dSel(IndexCoord(j))*pi/180);
        end
    end
    model_state.setQ(Qvect);
    model_state.setTime(dat.data(i,1));
    pause(dt);
    vis.show(model_state);
end


vis.delete();
model_state.delete();
clear('vis','model_state')


end