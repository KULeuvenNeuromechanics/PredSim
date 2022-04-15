% --------------------------------------------------------------------------
% test_muscleAnalysisAPI
%   This script tests the function "PreProcessing\muscleAnalysisAPI
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

%% compare results of API-based muscle analysis with analysis ran from xml file

[muscle_data_xml, Qs_dummy] = muscleAnalysis(S,osim_path,model_info,10);
muscle_spanning_joint_info_xml = squeeze(sum(abs(muscle_data_xml.dM), 1));
muscle_spanning_joint_info_xml(muscle_spanning_joint_info_xml<=0.0001 & muscle_spanning_joint_info_xml>=-0.0001) = 0;
muscle_spanning_joint_info_xml(muscle_spanning_joint_info_xml~=0) = 1;

muscle_spanning_joint_info_API = get_muscle_spanning_joint_info(S,osim_path,model_info);
model_info.muscle_info.muscle_spanning_joint_info = muscle_spanning_joint_info_API;
muscle_data_API = muscleAnalysisAPI(S,osim_path,model_info,Qs_dummy);


%%
% for i=1:length(muscle_data_xml.dof_names)
%     disp(muscle_data_xml.dof_names{i});
%     for j=1:length(muscle_data_xml.muscle_names)
%         if muscle_spanning_joint_info_xml(j,i)==1
%             disp(['   ' muscle_data_xml.muscle_names{j}])
%         end
%     end
% end
% 
% for i=1:length(muscle_data_API.dof_names)
%     disp(muscle_data_API.dof_names{i});
%     for j=1:length(muscle_data_API.muscle_names)
%         if muscle_spanning_joint_info_API(j,i)==1
%             disp(['   ' muscle_data_API.muscle_names{j}])
%         end
%     end
% end

%%
diff_msji = muscle_spanning_joint_info_xml - muscle_spanning_joint_info_API;
diff_lMT = muscle_data_xml.lMT - muscle_data_API.lMT;
diff_dM = muscle_data_xml.dM - muscle_data_API.dM;

diff_lMT_v = reshape(diff_lMT,size(diff_lMT,1)*size(diff_lMT,2),1);
diff_dM_v = reshape(diff_dM,size(diff_dM,1)*size(diff_dM,2)*size(diff_dM,3),1);

cond_lMT = norm(diff_lMT_v,"inf");
cond_dM = norm(diff_dM_v,"inf");







