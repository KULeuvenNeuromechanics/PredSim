
clear
close all
clc

[pathTests,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathTests);

addpath(fullfile(pathRepo,'OCP'))

load(fullfile(pathRepo,'Tests','ReferenceResults','Falisse_et_al_2022',['Falisse_et_al_2022','_paper.mat']),'R','model_info');


N = size(R.ground_reaction.GRF_r,1);

GRFk_opt0 = [R.ground_reaction.GRF_r,R.ground_reaction.GRF_l];

idx_1 = round(rand(1)*N);

GRFk_opt = GRFk_opt0([idx_1:N,1:idx_1-1],:);


[idx_GC,idx_GC_base_forward_offset,HS1,threshold] = getStancePhaseSimulation(GRFk_opt,model_info.mass/3);


GRFk_opt_GC = GRFk_opt(idx_GC,:);


figure
plot(GRFk_opt0)
hold on
plot(GRFk_opt_GC,'--')


