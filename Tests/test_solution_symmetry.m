
clear
close all
clc

[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
[pathRepoFolder,~,~] = fileparts(pathRepo);
addpath([pathRepo '/PostProcessing'])

%% General settings
% These settings will apply to all figures

% Construct a cell array with full paths to files with saved results for
% which you want to appear on the plotted figures.
results_folder = fullfile(pathRepoFolder,'PredSimResults');
subj = 'subject1_2D';
% subj = 'Fal_s1_mtp';
result_paths{1} = fullfile(results_folder,subj,[subj '_v6.mat']);

load(result_paths{1},'R','model_info')

%%

Qs1 = R.kinematics.Qs_rad(1:50,:);
Qs2 = R.kinematics.Qs_rad(51:100,:);
Qs2(:,model_info.ExtFunIO.symQs.QsInvA) = Qs2(:,model_info.ExtFunIO.symQs.QsInvB);
Qs2(:,model_info.ExtFunIO.symQs.QsOpp) = -Qs2(:,model_info.ExtFunIO.symQs.QsOpp);

diff_Qs = Qs2 - Qs1;

Qdots1 = R.kinematics.Qdots_rad(1:50,:);
Qdots2 = R.kinematics.Qdots_rad(51:100,:);
Qdots2(:,model_info.ExtFunIO.symQs.QdotsInvA) = Qdots2(:,model_info.ExtFunIO.symQs.QdotsInvB);
Qdots2(:,model_info.ExtFunIO.symQs.QsOpp) = -Qdots2(:,model_info.ExtFunIO.symQs.QsOpp);

diff_Qdots = Qdots2 - Qdots1;

Qddots1 = R.kinematics.Qddots_rad(1:50,:);
Qddots2 = R.kinematics.Qddots_rad(51:100,:);
Qddots2(:,model_info.ExtFunIO.symQs.QdotsInvA) = Qddots2(:,model_info.ExtFunIO.symQs.QdotsInvB);
Qddots2(:,model_info.ExtFunIO.symQs.QsOpp) = -Qddots2(:,model_info.ExtFunIO.symQs.QsOpp);

diff_Qddots = Qddots2 - Qddots1;

Ts1 = R.kinetics.T_ID(1:50,:);
Ts2 = R.kinetics.T_ID(51:100,:);
Ts2(:,model_info.ExtFunIO.symQs.QdotsInvA) = Ts2(:,model_info.ExtFunIO.symQs.QdotsInvB);
Ts2(:,model_info.ExtFunIO.symQs.QsOpp) = -Ts2(:,model_info.ExtFunIO.symQs.QsOpp);

diff_Ts = Ts2 - Ts1;

Ts01 = R.kinetics.T_ID_0(1:50,:);
Ts02 = R.kinetics.T_ID_0(51:100,:);
Ts02(:,model_info.ExtFunIO.symQs.QdotsInvA) = Ts02(:,model_info.ExtFunIO.symQs.QdotsInvB);
Ts02(:,model_info.ExtFunIO.symQs.QsOpp) = -Ts02(:,model_info.ExtFunIO.symQs.QsOpp);

diff_Ts0 = Ts02 - Ts01;

%%

import casadi.*
pathmain = pwd;
cd(R.S.misc.subject_path)
F  = external('F',R.S.misc.external_function);
cd(pathmain);

% [F1] = load_external_function(R.S);

N = size(Qs1,1);

QsQdots1 = zeros(N,2*model_info.ExtFunIO.jointi.nq.all);
QsQdots1(:,1:2:end) = Qs1;
QsQdots1(:,2:2:end) = Qdots1;

QsQdots2 = zeros(N,2*model_info.ExtFunIO.jointi.nq.all);
QsQdots2(:,1:2:end) = Qs2;
QsQdots2(:,2:2:end) = Qdots2;

Foutk_opt1 = zeros(N,F.nnz_out);
Foutk_opt2 = Foutk_opt1;

for i = 1:N
    % ID moments
    [res1] = F([QsQdots1(i,:)';Qddots1(i,:)']);
    Foutk_opt1(i,:) = full(res1);

    [res2] = F([QsQdots2(i,:)';Qddots2(i,:)']);
    Foutk_opt2(i,:) = full(res1);
end

% Foutk_opt1 = Foutk_opt(:,1:model_info.ExtFunIO.jointi.nq.all);
% Foutk_opt2 = Foutk_opt1(idx_gc_inv,:);
Foutk_opt2(:,model_info.ExtFunIO.symQs.QdotsInvA) = Foutk_opt2(:,model_info.ExtFunIO.symQs.QdotsInvB);
Foutk_opt2(:,model_info.ExtFunIO.symQs.QsOpp) = -Foutk_opt2(:,model_info.ExtFunIO.symQs.QsOpp);
% Foutk_opt2 = Foutk_opt2(idx_gc_inv,:);

diff = Foutk_opt2 - Foutk_opt1;

