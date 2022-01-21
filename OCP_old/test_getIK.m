clear
clc

% add repository to workspace
% [pathHere,~,~] = fileparts(mfilename('fullpath'));
% [pathRepo,~,~] = fileparts(pathHere);
% addpath([pathRepo '/PreProcessing'])
% % addpath([pathRepo '/PreProcessing'])

% % load model info
% load(fullfile(pathRepo,'Subjects/test_1/default_model_info.mat'))

osim_path = 'test_1.osim';
% % load model info
% load('model_info.mat')

%% make tmp model_info
load('F_test_1_IO.mat');

model_info.ExtFunIO = IO;

import org.opensim.modeling.*;
model = Model(osim_path);
for i=1:model.getMuscles().getSize()
    muscle_names{i} = char(model.getMuscles().get(i-1).getName());
end

model_info.muscle_info.muscle_names = muscle_names;

%% create motion file
% colheaders in random order
coordinate_names = fieldnames(model_info.ExtFunIO.coordi);
coordinate_names_random = coordinate_names(randperm(length(coordinate_names)));
q_in.labels = [{'time'},coordinate_names_random'];
% random data
q_in.data = [linspace(0.01,0.05,10)', -pi/180*rand([10,length(coordinate_names)])];
% angles are in radians
q_in.inDeg = 'no';
% filename
filename = 'test_IK.mot';
% generate motion file
write_motionFile_v40(q_in, filename)

%%
joints = fieldnames(model_info.ExtFunIO.coordi)';
Qs = getIK(filename,model_info);
Qs_2 = getIK_original(filename,joints);

%% compare getIK to data

for i=1:length(coordinate_names)
    coord_name_i = coordinate_names{i};
    res = q_in.data(:,strcmp(q_in.labels,coord_name_i)) - Qs.all(:,strcmp(Qs.colheaders,coord_name_i));
    disp([coord_name_i '    ' num2str(max(abs(res)))])
end

%% compare getIK to getIK original

for i=1:length(Qs.colheaders)
    coord_name_i = Qs.colheaders{i};
    res = Qs_2.all(:,i) - Qs.all(:,i);
    disp([coord_name_i '    ' num2str(max(abs(res)))])
end
