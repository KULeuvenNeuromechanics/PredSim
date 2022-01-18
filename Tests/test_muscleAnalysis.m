%
% Verification tests of muscleAnalysis.m
%
% Author: Lars D'Hondt
% Date: 17/January/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

%% preparation
% paths
[path_tests,~,~] = fileparts(mfilename('fullpath'));
[path_main,~,~] = fileparts(path_tests);
addpath([path_main '\PreProcessing'])

% settings
S.path_main = path_main;
S.subject.name = 'test_1';
osim_path = fullfile(S.path_main,'Subjects',S.subject.name,[S.subject.name '.osim']);

% model_info
load(fullfile(S.path_main,'Subjects\test_1\F_test_1_IO.mat'));

model_info.ExtFunIO = IO;

import org.opensim.modeling.*;
model = Model(osim_path);
for i=1:model.getMuscles().getSize()
    muscle_names{i} = char(model.getMuscles().get(i-1).getName());
end

model_info.muscle_info.muscle_names = muscle_names;

%% run function
% takes 12-ish minutes

% muscleAnalysis
% [MuscleData] = muscleAnalysis(S,osim_path,model_info);

% PolynomialFit
% [muscle_spanning_joint_info,MuscleInfo] = PolynomialFit(MuscleData);

%% save results
% save(fullfile(S.path_main,'Subjects',S.subject.name,'MuscleData.mat'),'MuscleData')
% save(fullfile(S.path_main,'Subjects',S.subject.name,'muscle_spanning_joint_info.mat'),'muscle_spanning_joint_info')
% save(fullfile(S.path_main,'Subjects',S.subject.name,'MuscleInfo.mat'),'MuscleInfo');

%% load results
% reference
% load(fullfile(path_tests,'reference_muscleAnalysis','MuscleData.mat'),'MuscleData')
load(fullfile(path_tests,'reference_muscleAnalysis','muscle_spanning_joint_INFO.mat'),'muscle_spanning_joint_INFO')
load(fullfile(path_tests,'reference_muscleAnalysis','MuscleInfo.mat'),'MuscleInfo');
muscle_spanning_joint_info_ref = muscle_spanning_joint_INFO;
MuscleInfo_ref = MuscleInfo;

% new 
% load(fullfile(S.path_main,'Subjects',S.subject.name,'MuscleData.mat'),'MuscleData')
load(fullfile(S.path_main,'Subjects',S.subject.name,'muscle_spanning_joint_info.mat'),'muscle_spanning_joint_info')
load(fullfile(S.path_main,'Subjects',S.subject.name,'MuscleInfo.mat'),'MuscleInfo');

%% compare right side to reference

for i=1:length(MuscleInfo_ref.muscle)
    idx = find(strcmp(MuscleInfo_ref.muscle(i).m_name,{MuscleInfo.muscle(:).m_name}));
    test_dofs(i) = min(strcmp(MuscleInfo_ref.muscle(i).DOF,MuscleInfo.muscle(idx).DOF'));
    test_coeffs(i) = max(abs( (MuscleInfo_ref.muscle(i).coeff - MuscleInfo.muscle(idx).coeff)) );
end


%% check symmetry

% MuscleInfo_right = cell(0);
% MuscleInfo_left = cell(0);
n_r = 1;
n_l = 1;

for i=1:length(MuscleInfo.muscle)
    m_name_i = MuscleInfo.muscle(i).m_name;
    if strcmp(m_name_i(end),'r')
        MuscleInfo_right(n_r) = MuscleInfo.muscle(i);
        n_r = n_r+1;
    elseif strcmp(m_name_i(end),'l')
        MuscleInfo_left(n_l) = MuscleInfo.muscle(i);
        n_l = n_l+1;
    end
end

for i=1:length(MuscleInfo_right)
    m_name_i = MuscleInfo_right(i).m_name;
    m_name_i(end) = 'l';
    idx = find(strcmp(m_name_i,{MuscleInfo_left(:).m_name}));

    
    test_coeffs_symm(i) = max(abs( (MuscleInfo_right(i).coeff - MuscleInfo_left(idx).coeff)) );

    for j=1:length(MuscleInfo_right(i).DOF)
        dof_j = MuscleInfo_right(i).DOF{j};
        if strcmp(dof_j(end),'r')
            dof_j(end) = 'l';
        end

        test_dofs_symm_i(j) = strcmp(dof_j,MuscleInfo_left(idx).DOF{j});

    end
    test_dofs_symm(i) = min(test_dofs_symm_i);


end






