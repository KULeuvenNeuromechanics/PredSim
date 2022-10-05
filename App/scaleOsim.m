function scaleOsim(MainPath, U)

MainPath = 'C:\GBW_MyPrograms\PredSim_test';
U.ModelName = 'test_scale';
U.Mass = 75;
U.Height = 85;

%% inputs

osim_output_name = [U.ModelName '.osim'];
output_dir = fullfile(MainPath, 'Subjects',U.ModelName);
if ~isfolder(output_dir)
    mkdir(output_dir);
end

generic_scale = fullfile(MainPath, 'OpenSimModel', 'uniformScale.xml');


%% adjust scaling setup file

scalefile = xmlread(generic_scale);

% change model name
scalefile.getElementsByTagName('ScaleTool').item(0).getAttributes.item(0).setValue(U.ModelName)

% change mass
scalefile.getElementsByTagName('mass').item(0).getFirstChild.setNodeValue(num2str(U.Mass));

% change scaling factors
scaleF = U.Height/1.80;
segementscaling = [num2str(scaleF) ' ' num2str(scaleF) ' ' num2str(scaleF)];

%get amount of segments
nSegments = scalefile.getElementsByTagName('Scale').getLength;

for i=1:nSegments
    %adjust scales
    
    scalefile.getElementsByTagName('Scale').item(i-1).getElementsByTagName('scales').item(0).getFirstChild.setNodeValue(segementscaling);
end

% set output_model_file

scalefile.getElementsByTagName('output_model_file').item(0).getFirstChild.setNodeValue(osim_output_name)

%% export adjusted scaling setup file
output_name = fullfile(output_dir,[U.ModelName '_Scale_Setup.xml']);
xmlwrite(output_name,scalefile);

%% scale generic model using adjusted scaling setup file

import org.opensim.modeling.*

scale = ScaleTool(output_name)
scale.run()

end

