% Script to test the function getIK
%
% Author: Lars D'Hondt
% Date: 1 dec 2021
%
%--------------------------------------------------------------------------
clear
clc

% add repository to workspace
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/PreProcessing'])
% addpath([pathRepo '/PreProcessing'])

% load model info
load(fullfile(pathRepo,'Subjects/test_1/default_model_info.mat'))

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
filename = fullfile(pathHere,'test_IK.mot');
% generate motion file
write_motionFile_v40(q_in, filename)

%%
Qs = getIK(filename,model_info);

%%


for i=1:length(coordinate_names)
    coord_name_i = coordinate_names{i};
    res = q_in.data(:,strcmp(q_in.labels,coord_name_i)) - Qs.all(:,strcmp(Qs.colheaders,coord_name_i));
    disp([coord_name_i '    ' num2str(max(abs(res)))])
end













