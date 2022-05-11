
clear
close all
clc


%% reference result
ResultsRepo = 'C:\Users\u0150099\OneDrive - KU Leuven\3dpredictsim_results';
ref_file = fullfile([ResultsRepo '/debug\Fal_s1_mtp_sd_MTPp_k17_d05_ig24_pp.mat']);
load(ref_file,'R')
R_ref = R;
x = 1:(100-1)/(size(R_ref.Qs,1)-1):100;


load('C:\Users\u0150099\Documents\PredSimResults\Fal_s1_mtp\Fal_s1_mtp_v1.mat','R','model_info')
R_mtp = R;
model_info_mtp = model_info;

load('C:\Users\u0150099\Documents\PredSimResults\Fal_s1\Fal_s1_v2.mat','R','model_info')

figure
for i=1:31
    subplot(4,8,i)
    hold on
    plot(x,R_ref.Qs(:,i))
    title(R_ref.colheaders.joints{i},'Interpreter','none')
end

for i=1:31
    idx = model_info.ExtFunIO.coordi.(R_ref.colheaders.joints{i});
    subplot(4,8,idx)
    hold on
    plot(x,R.Qs(:,i))
    plot(x,R_mtp.Qs(:,i))
end
legend('reference','new code','new code mtp')

%%
figure
idx_r = [];
for i=1:length(R_ref.colheaders.muscles)
    if strcmp(R_ref.colheaders.muscles{i}(end),'r')
        idx_r(end+1) = i;
    end
end
for i=1:length(idx_r)
    subplot(7,7,i)
    hold on
    plot(x,R_ref.a(:,idx_r(i)))
    title(R_ref.colheaders.muscles{idx_r(i)},'Interpreter','none')
end

idx_r2 = [];
for i=1:length(model_info.muscle_info.muscle_names)
    if strcmp(model_info.muscle_info.muscle_names{i}(end),'r')
        idx_r2(end+1) = i;
    end
end
for i=1:length(idx_r2)
    idx0 = find(strcmp(model_info.muscle_info.muscle_names(idx_r2),R_ref.colheaders.muscles(idx_r(i))));
%     idx = find(idx_r(:)==idx0);
    subplot(7,7,i)
    hold on
    plot(x,R.a(:,idx_r2(idx0)))
    plot(x,R_mtp.a(:,idx_r2(idx0)))
end
legend('reference','new code','new model mtp')


%%

% clc
% for i=1:length(model_info.ExtFunIO.symQs.QsInvA)
%     name1 = model_info.ExtFunIO.coord_names.all{model_info.ExtFunIO.symQs.QsInvA(i)};
%     name2 = model_info.ExtFunIO.coord_names.all{model_info.ExtFunIO.symQs.QsInvB(i)};
%     disp([name1 '  -  ' name2])
% end


