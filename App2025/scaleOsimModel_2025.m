function scaleOsimModel_2025(MainPath, U)

% % test inputs
% MainPath = 'C:\GBW_MyPrograms\PredSim_test';
% U.ModelName = 'test_scale';
% U.Height = 1.80;
% 
% sf.foot = 1;
% sf.low_leg = 1;
% sf.upp_leg = 1;
% sf.upp_arm = 1;
% sf.low_arm = 1;
% sf.shoulder = 1;
% 
% % end of test inputs




%% inputs
% folders
osim_output_name = [U.ModelName '.osim'];
output_dir = fullfile(MainPath, 'Subjects',U.GroupName, U.ModelName);
if exist(fullfile(output_dir,osim_output_name),'file')
    disp('Scaled model already exists.')
    return
end

if ~isfolder(output_dir)
    mkdir(output_dir);
end

% copy unscaled model file
copyfile(fullfile(MainPath,'Subjects','gait1018','gait1018.osim'),output_dir);

% generic xml file
generic_scalefile = fullfile(MainPath,'Subjects','Vitruvian_Man','uniformScale.xml');

% general scale factor
sf_total = U.Height/160;

% apply scale factors to segments
sf_seg.pelvis = [1,1,1]*sf_total;
sf_seg.femur_r = [1,1,1]*sf_total;
sf_seg.femur_l = [1,1,1]*sf_total;
sf_seg.tibia_r = [1,1,1]*sf_total;
sf_seg.tibia_l = [1,1,1]*sf_total;
sf_seg.talus_r = [1,1,1]*sf_total;
sf_seg.calcn_r = [1,1,1]*sf_total;
sf_seg.toes_r = [1,1,1]*sf_total;
sf_seg.talus_l = [1,1,1]*sf_total;
sf_seg.calcn_l = [1,1,1]*sf_total;
sf_seg.toes_l = [1,1,1]*sf_total;
sf_seg.torso = [sqrt(sf_total*sf_total),sf_total,sf_total];
sf_seg.humerus_r = [1,1,1]*sf_total;
sf_seg.humerus_l = [1,1,1]*sf_total;
sf_seg.ulna_r = [1,1,1]*sf_total;
sf_seg.radius_r = [1,1,1]*sf_total;
sf_seg.hand_r = [1,1,1]*sf_total;
sf_seg.ulna_l = [1,1,1]*sf_total;
sf_seg.radius_l = [1,1,1]*sf_total;
sf_seg.hand_l = [1,1,1]*sf_total;

%% adjust scaling setup file

scalefile = xmlread(generic_scalefile);

% change model name
scalefile.getElementsByTagName('ScaleTool').item(0).getAttributes.item(0).setValue(U.ModelName)

% change mass
scalefile.getElementsByTagName('mass').item(0).getFirstChild.setNodeValue(num2str(U.Mass));

% set model_file
scalefile.getElementsByTagName('model_file').item(0).getFirstChild.setNodeValue('gait1018.osim');

% get amount of segments
nSegments = scalefile.getElementsByTagName('Scale').getLength;

for i=1:nSegments
    % get name of segment
    segment_name = scalefile.getElementsByTagName('Scale').item(i-1).getElementsByTagName('segment').item(0).getFirstChild.getNodeValue;

    % segment-specific scale factor
    if isfield(sf_seg,char(segment_name))
        sf_segm = sf_seg.(char(segment_name));
    else
        sf_segm = [1,1,1];
    end
    % 
    % % add general scale factor
    % segm_sf = segm_sf*sf_total;
    % %
    segment_scaling = [num2str(sf_segm(1)) ' ' num2str(sf_segm(2)) ' ' num2str(sf_segm(3))];

    %adjust scales
    scalefile.getElementsByTagName('Scale').item(i-1).getElementsByTagName('scales').item(0).getFirstChild.setNodeValue(segment_scaling);
end

% set output_model_file
scalefile.getElementsByTagName('output_model_file').item(0).getFirstChild.setNodeValue(osim_output_name)

%% export adjusted scaling setup file
output_name = fullfile(output_dir,[U.ModelName '_Scale_Setup.xml']);
xmlwrite(output_name,scalefile);

%% scale generic model using adjusted scaling setup file

import org.opensim.modeling.*

scale = ScaleTool(output_name);
scale.run();


%% adapt contact sphere placement
fixContactSpherePositionAfterScaling(fullfile(output_dir,'gait1018.osim'),fullfile(output_dir,osim_output_name));
scaleContactSpheres_2025(fullfile(output_dir,'gait1018.osim'),fullfile(output_dir,osim_output_name));

% remove unscaled model file
delete(fullfile(output_dir,'gait1018.osim'));

end