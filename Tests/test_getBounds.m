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
addpath([pathRepo '/OCP'])
addpath([pathRepo '/VariousFunctions'])


%% create needed settings
S.subject.folder_name = pathHere;
S.subject.IK_Bounds = 'test_IK.mot';
S.subject.v_pelvis_x_trgt = 1.25;
S.misc.gaitmotion_type = 'HalfGaitCycle';
S.bounds.coordinates = [];
S.bounds.a.lower = 0.05;
S.bounds.t_final.lower = 0.5;
S.bounds.t_final.upper = 2;

%% create model_info
% some fields are required to have a value
load(fullfile(pathRepo,'Subjects/Hamner_modified/F_Hamner_modified_IO.mat'))
model_info.ExtFunIO = IO;
model_info.IG_pelvis_y = 0.87;
model_info.muscle_info.NMuscle = 30;
model_info.muscle_info.tdeact = 0.5;
model_info.muscle_info.tact = 0.5;
model_info.ExtFunIO.jointi.nq.torqAct = 0;
model_info.actuator_info.parameters = [];

coordinate_names = fieldnames(model_info.ExtFunIO.coordi);
model_info.ExtFunIO.coord_names.all = coordinate_names;
model_info.ExtFunIO.jointi.nq.all = length(coordinate_names);
model_info.ExtFunIO.jointi.base_forward = 4;
model_info.ExtFunIO.jointi.base_lateral = 6;
tmp = 1:length(coordinate_names);
model_info.ExtFunIO.symQs.QsInvA = tmp;
model_info.ExtFunIO.symQs.QsInvB = tmp;
model_info.ExtFunIO.symQs.QdotsInvA = tmp;
model_info.ExtFunIO.symQs.QdotsInvB = tmp;

%% create motion file
% colheaders
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









