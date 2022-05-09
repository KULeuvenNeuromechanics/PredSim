
clear
close all
clc


%% reference result
ResultsRepo = 'C:\Users\u0150099\OneDrive - KU Leuven\3dpredictsim_results';
ref_file = fullfile([ResultsRepo '/debug\Fal_s1_mtp_sd_MTPp_k17_d05_ig24_pp.mat']);
load(ref_file,'R')
R_ref = R;
x = 1:(100-1)/(size(R_ref.Qs,1)-1):100;

load('C:\Users\u0150099\Documents\PredSimResults\Fal_s1\R_Fal_s1.mat','R','model_info')

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
end
legend('reference','new code')
