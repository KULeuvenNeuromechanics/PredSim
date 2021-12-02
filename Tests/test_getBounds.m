% Script to test the function getBounds
%
% Author: Lars D'Hondt
% Date: 2 dec 2021
%
% TO DO
%--------------------------------------------------------------------------
clear
clc

% add repository to workspace
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/PreProcessing'])

% load model info
load(fullfile(pathRepo,'Subjects/test_1/default_model_info.mat'))

%% create motion file
% colheaders
coordinate_names = fieldnames(model_info.ExtFunIO.coordi);
coordinate_names_random = coordinate_names(randperm(length(coordinate_names)));
q_in.labels = [{'time'},coordinate_names'];
% random data
q_in.data = [linspace(0.01,0.05,10)', -pi/180*(rand([10,length(coordinate_names)])-0.3)*0.5];
% angles are in radians
q_in.inDeg = 'no';
% filename
filename = fullfile(pathHere,'test_IK.mot');
% generate motion file
write_motionFile_v40(q_in, filename)

%% create needed settings
S.subject.folder_name = pathHere;
S.subject.IKfile_bounds = 'test_IK.mot';
S.subject.vPelvis_x_trgt = 1.25;
S.subject.IG_PelvisY = 0.85;
S.subject.mtp_type = 'active';

%% call getBounds
[bounds,scaling] = getBounds(S,model_info);

%%
for i=1:length(coordinate_names)
    coord_name_i = coordinate_names{i};
    q = q_in.data(:,strcmp(q_in.labels,coord_name_i));
    lb = bounds.Qs.lower(i);
    ub = bounds.Qs.upper(i);
    sf = scaling.Qs(i);
    [~,idx] = max(abs(q));
    res = q(idx)/sf;
    if lb<res && res<ub
        check = 'ok';
    else
        check = 'nok';
    end

    disp([coord_name_i '    ' num2str(lb) ' < ' num2str(res) ' < ' num2str(ub) '  ' check])
end









