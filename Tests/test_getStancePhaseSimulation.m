
clear
close all
clc

addpath('../OCP')

load('.\Falisse_et_al_2022_Results\Falisse_et_al_2022_v1.mat','R','model_info');

N = size(R.ground_reaction.GRF_r,1);

GRFk_opt0 = [R.ground_reaction.GRF_r,R.ground_reaction.GRF_l];

idx_1 = round(rand(1)*N);

GRFk_opt = GRFk_opt0([idx_1:N,1:idx_1-1],:);


[idx_GC,HS1,threshold] = getStancePhaseSimulation(GRFk_opt,N,model_info.mass/3);


GRFk_opt_GC = GRFk_opt(idx_GC,:);


figure
plot(GRFk_opt0)
hold on
plot(GRFk_opt_GC,'--')


