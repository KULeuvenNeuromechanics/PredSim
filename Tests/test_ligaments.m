
clear
clc

import casadi.*


load('C:\GBW_MyPrograms\PredSimResults\DHondt_2023_3seg\DHondt_2023_3seg_job666.mat')


f_mtjLigament = Function.load(fullfile('C:\Users\u0150099\Documents\master_thesis\3dpredictsim\Polynomials\Fal_s1_mtjc4_FK_sc_FDB2_MTc5','f_getMtjLigamentMoment_exp5'));

%%

Qs = MX.sym('Qs',model_info.ExtFunIO.jointi.nq.all,1);
figure
tiledlayout('flow')
nexttile
spy(jacobian(f_casadi.ligamentMoment(Qs),Qs));

nexttile
spy(jacobian(f_casadi.ligamentMoment_single(Qs),Qs));

nexttile
spy(jacobian(f_casadi.ligamentMoment_multi(Qs),Qs));

%%

import casadi.*

Qs = MX(model_info.ExtFunIO.jointi.nq.all,1);
q_mtj = MX.sym('q_mtj',1,1);
Qs(model_info.ExtFunIO.coordi.mtj_angle_r,1) = q_mtj;

M_lig = f_casadi.ligamentMoment_single(Qs);
M_mtj = M_lig(model_info.ExtFunIO.coordi.mtj_angle_r);

f_mtj_lig = Function('f_mtj_lig',{q_mtj},{M_mtj});

% diff = compareCasADiFunctions(f_mtjLigament,f_mtj_lig,[],10,[-3,0])


qs = linspace(-20,20,200);
M1 = full(f_mtjLigament(qs*pi/180));
M2 = full(f_mtj_lig(qs*pi/180));

figure
hold on
plot(qs,M1)
plot(qs,M2)




%%

ResultsRepo = 'C:\Users\u0150099\OneDrive - KU Leuven\3dpredictsim_results';
ref = load(fullfile([ResultsRepo '\with_better_knee\Fal_s1_mtjc4_FK_sc_cspx10_cg9_o1x10_ATx50_TFMox120_Fpsl10_MTc5_MTPm_k1_d01_tau_MTJm_nl_MG_exp5_table_d01_PF_Natali2010_ls146_FDB2_lTs125_Fpsl10_ig21.mat']),'setup');

%%
% Qs_ub = setup.bounds.Qs.upper.*setup.scaling.Qs;
% Qs_lb = setup.bounds.Qs.lower.*setup.scaling.Qs;
% 
% Qs_ub_ref = ref.setup.bounds.Qs.upper.*ref.setup.scaling.Qs;
% Qs_lb_ref = ref.setup.bounds.Qs.lower.*ref.setup.scaling.Qs;
% 
% 
% find(abs([Qs_ub_ref-Qs_ub]')>1e-6)
% 
% find(abs([Qs_lb_ref-Qs_lb]')>1e-6)

%%

ref = load(fullfile([ResultsRepo '\results_paper\Fal_s1_mtjc4_FK_sc_cspx10_cg9_o1x10_ATx50_TFMox120_Fpsl10_MTc5_MTPm_k1_d01_tau_MTJm_nl_MG_exp5_table_d01_PF_Natali2010_ls146_FDB2_lTs125_Fpsl10_ig1_N100_pp.mat']));


Qs_ref = ref.R.Qs;
M_lig_ref = zeros(size(Qs_ref));

for i=1:size(Qs_ref,1)
    Qs_i = Qs_ref(i,:)*pi/180;
    M_i = f_casadi.ligamentMoment_multi(Qs_i);

    M_lig_ref(i,:) = full(M_i);
    
end


M_mtj_PF_ref = ref.R.windlass.MA_PF.mtj.*ref.R.windlass.F_PF;
M_mtp_PF_ref = ref.R.windlass.MA_PF.mtp.*ref.R.windlass.F_PF;

figure
tiledlayout('flow')
nexttile
hold on
plot(M_mtj_PF_ref)
plot(-M_lig_ref(:,model_info.ExtFunIO.coordi.mtj_angle_r))

nexttile
hold on
plot(M_mtp_PF_ref)
plot(-M_lig_ref(:,model_info.ExtFunIO.coordi.mtp_angle_r))


%%

CSA = 60;
ls = 0.146;
l = linspace(0.98,1.1,200)*ls;
F = plantarFasciaNatali2010(CSA,ls,l);

figure
hold on
plot(l,F)

