
clear
clc

% paths
[pathTests,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathTests);
addpath([pathRepo '\PreProcessing'])
addpath([pathRepo '\VariousFunctions'])

load(fullfile(pathTests,'Falisse_et_al_2022_Results\Falisse_et_al_2022_v1.mat'))

S = R.S;
osim_path = replace(model_info.osim_path, S.misc.main_path, pathRepo);
% osim_path = model_info.osim_path;

S.subject.base_joints_legs = 'hip';
S.subject.base_joints_arms = {'acromial_r'};

[symQs] = getCoordinateSymmetry(S,osim_path,model_info);

[symmetry, jointi] = identify_kinematic_chains(S,osim_path,model_info);


max(abs(symQs.QsInvA - symmetry.QsInvA))
max(abs(symQs.QsInvB - symmetry.QsInvB))
max(abs(symQs.QdotsInvA - symmetry.QdotsInvA))
max(abs(symQs.QdotsInvB - symmetry.QdotsInvB))



