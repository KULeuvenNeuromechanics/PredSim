%
% Author: Lars D'Hondt
% Date: 17/January/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

[path_tests,~,~] = fileparts(mfilename('fullpath'));
[path_main,~,~] = fileparts(path_tests);
addpath([path_main '\PreProcessing'])
S.path_main = path_main;
S.subject.name = 'test_1';
osim_path = fullfile(S.path_main,'Subjects',S.subject.name,[S.subject.name '.osim']);

%% make tmp model_info
load(fullfile(S.path_main,'Subjects\test_1\F_test_1_IO.mat'));

model_info.ExtFunIO = IO;

import org.opensim.modeling.*;
model = Model(osim_path);
for i=1:model.getMuscles().getSize()
    muscle_names{i} = char(model.getMuscles().get(i-1).getName());
end

model_info.muscle_info.muscle_names = muscle_names;

%%
% [MuscleData] = muscleAnalysis(S,osim_path,model_info);

%% Call PolynomialFit
% [muscle_spanning_joint_info,MuscleInfo] = PolynomialFit(MuscleData);

% save(fullfile(S.path_main,'Subjects',S.subject.name,'MuscleData.mat'),'MuscleData')
% save(fullfile(S.path_main,'Subjects',S.subject.name,'muscle_spanning_joint_info.mat'),'muscle_spanning_joint_info')
% save(fullfile(S.path_main,'Subjects',S.subject.name,'MuscleInfo.mat'),'MuscleInfo');

%%
% load(fullfile(path_tests,'reference_muscleAnalysis','MuscleData.mat'),'MuscleData')
load(fullfile(path_tests,'reference_muscleAnalysis','muscle_spanning_joint_INFO.mat'),'muscle_spanning_joint_INFO')
load(fullfile(path_tests,'reference_muscleAnalysis','MuscleInfo.mat'),'MuscleInfo');
muscle_spanning_joint_info_ref = muscle_spanning_joint_INFO;
MuscleInfo_ref = MuscleInfo;

% load(fullfile(S.path_main,'Subjects',S.subject.name,'MuscleData.mat'),'MuscleData')
load(fullfile(S.path_main,'Subjects',S.subject.name,'muscle_spanning_joint_info.mat'),'muscle_spanning_joint_info')
load(fullfile(S.path_main,'Subjects',S.subject.name,'MuscleInfo.mat'),'MuscleInfo');

%%

for i=1:length(MuscleInfo_ref.muscle)
    idx = find(strcmp(MuscleInfo_ref.muscle(i).m_name,{MuscleInfo.muscle(:).m_name}));
    test_dofs(i) = min(strcmp(MuscleInfo_ref.muscle(i).DOF,MuscleInfo.muscle(idx).DOF'));
    test_coeffs(i) = max(abs( (MuscleInfo_ref.muscle(i).coeff - MuscleInfo.muscle(idx).coeff)) );

end

