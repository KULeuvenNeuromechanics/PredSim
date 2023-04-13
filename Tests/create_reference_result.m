clear
clc


ResultsRepo = 'C:\Users\u0150099\OneDrive - KU Leuven\3dpredictsim_results';

% load(fullfile([ResultsRepo '\results_paper\Fal_s1_mtp_FK_sc_cspx10_cg9_o1x10_ATx50_TFMox120_Fpsl10_MTc5_MTPp_k25_d020_tau_ig1_N100_pp.mat']),'R');

load(fullfile([ResultsRepo '\results_paper\Fal_s1_mtjc4_FK_sc_cspx10_cg9_o1x10_ATx50_TFMox120_Fpsl10_MTc5_MTPm_k1_d01_tau_MTJm_nl_MG_exp5_table_d01_PF_Natali2010_ls146_FDB2_lTs125_Fpsl10_ig1_N100_pp.mat']),'R');

R_ref = R;
clearvars('R')

R.kinematics.Qs = R_ref.Qs;
R.kinetics.T_ID = R_ref.Tid;
R.muscles.a = R_ref.a;
R.colheaders.coordinates = R_ref.colheaders.joints;
R.colheaders.muscles = R_ref.colheaders.muscles;

R.S.post_process.result_filename = 'ref_3segment';

N = size(R_ref.Qs,1);
NMuscle = size(R_ref.a,2);
R.S.metabolicE.model = 'Bhargava2004';
R.metabolics.Bhargava2004.Edot_gait = R_ref.MetabB.Etot;
R.metabolics.Bhargava2004.Adot = R_ref.MetabB.Adot;
R.metabolics.Bhargava2004.Mdot = R_ref.MetabB.Mdot;
R.metabolics.Bhargava2004.Sdot = R_ref.MetabB.Sdot;
R.metabolics.Bhargava2004.Wdot = R_ref.MetabB.Wdot;
R.metabolics.Bhargava2004.Edot_incl_basal = zeros(N,1);

R.muscles.FT = R_ref.FT;
R.muscles.Fce = R_ref.Muscle.Fce;
R.muscles.Fpass = R_ref.Muscle.Fpas;
% % R.muscles.Fiso = 
% R.muscles.lM = 
R.muscles.lMtilde = R_ref.lMtilde;
R.muscles.vM = R_ref.Muscle.vM;
R.muscles.vMtilde = R_ref.vMtilde;
R.muscles.lT = R_ref.lT;
R.muscles.vT = R_ref.vT;


model_info = [];

% R2path = 'C:\GBW_MyPrograms\PredSimResults\DHondt_2023_2seg\DHondt_2023_2seg_ref.mat';
R2path = 'C:\GBW_MyPrograms\PredSimResults\DHondt_2023_3seg\DHondt_2023_3seg_ref.mat';

save(R2path,'R','model_info')



