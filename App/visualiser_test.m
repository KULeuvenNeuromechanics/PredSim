% initialise the opensim visualisation
import org.opensim.modeling.*
[pathApp,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathApp);
modelPath = fullfile(pathRepo,'Subjects','Vitruvian_Man','Vitruvian_Man.osim');
osimModel = Model(modelPath);
osimModel.setUseVisualizer(true);
model_state = osimModel.initSystem();
vis = osimModel.getVisualizer();
% vis.addDirToGeometrySearchPaths(char(fullfile(S.pathRepo,'App','Geometry')));
vis.addDirToGeometrySearchPaths('C:\GBW_MyPrograms\OpenSim 4.3\Geometry');
pause(0.1);
vis.show(model_state);


