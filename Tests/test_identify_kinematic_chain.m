
clear
clc

load('C:\GBW_MyPrograms\PredSim_test\Tests\Falisse_et_al_2022_Results\Falisse_et_al_2022_v1.mat')

S = R.S;
osim_path = model_info.osim_path;

S.subject.base_joints_legs = 'hip';
S.subject.base_joints_arms = {'acromial_r'};

[symQs] = getCoordinateSymmetry(S,osim_path,model_info);

[symmetry, jointi] = identify_kinematic_chains(S,osim_path,model_info);


max(abs(symQs.QsInvA - symmetry.QsInvA))
max(abs(symQs.QsInvB - symmetry.QsInvB))
max(abs(symQs.QdotsInvA - symmetry.QdotsInvA))
max(abs(symQs.QdotsInvB - symmetry.QdotsInvB))

