% --------------------------------------------------------------------------
% 
% Original author: Lars D'Hondt
% Original date: 12/April/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

clear
close all
clc

%% preparation
% paths
[pathTests,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathTests);
addpath([pathRepo '\PreProcessing'])

S.misc.main_path = pathRepo;

S.subject.name = 'test_1'; 
S.subject.save_folder = fullfile(pathRepo,'test_1');
osim_path = fullfile(pathRepo,'Subjects','test_1','test_1.osim');
S.misc.subject_path = fullfile(S.misc.main_path,'Subjects',S.subject.name);
model_info = get_model_info(S,osim_path);

%%
muscle_spanning_joint_info = get_muscle_spanning_joint_info(S,osim_path,model_info);
model_info.muscle_info.muscle_spanning_joint_info = muscle_spanning_joint_info;
% muscle_data = muscleAnalysisAPI(S,osim_path,model_info,Qs_dummy);



msji1 = muscle_spanning_joint_info(model_info.ExtFunIO.symQs.orderMus,model_info.ExtFunIO.symQs.QsInvA);

msji2 = muscle_spanning_joint_info(model_info.ExtFunIO.symQs.orderMusInv,model_info.ExtFunIO.symQs.QsInvB);

diff_msji = norm(msji1-msji2,"inf")

