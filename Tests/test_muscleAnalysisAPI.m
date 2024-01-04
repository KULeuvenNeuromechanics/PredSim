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
addpath([pathRepo '\VariousFunctions'])
addpath([pathRepo '\DefaultSettings'])

load(fullfile(pathRepo,'Tests','ReferenceResults','Falisse_et_al_2022',['Falisse_et_al_2022','_paper.mat']),'R','model_info');
S = R.S;
osim_path = replace(model_info.osim_path, S.misc.main_path, pathRepo);
model_info.osim_path = osim_path;
S.subject.IG_selection = replace(S.subject.IG_selection, S.misc.main_path, pathRepo);
S.subject.IK_Bounds = replace(S.subject.IK_Bounds, S.misc.main_path, pathRepo);
S.misc.main_path = pathRepo;
S.misc.msk_geom_n_samples = 10;

[S] = getDefaultSettings(S,osim_path);

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

cond_lMT = max(abs(diff_lMT_v),[],'all');
cond_dM = max(abs(diff_dM_v),[],'all');

if cond_lMT > S.misc.threshold_lMT_fit
    fprintf('muscle-tendon lengths not ok\n')
else
    fprintf('muscle-tendon lengths ok\n')
end
if cond_dM > S.misc.threshold_dM_fit
    fprintf('muscle moment arms not ok\n')
else
    fprintf('muscle moment arms ok\n')
end




