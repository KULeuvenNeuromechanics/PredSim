
clear
close all
clc


%% reference result
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
ResultsRepo = [pathRepo '\test_results'];
addpath([pathRepo '\PlotFigures'])
ref_file = fullfile([ResultsRepo '\Fal_s1_mtp\Fal_s1_mtp_v14.mat']);
load(ref_file,'R','model_info');


coord_muscle_names_sel = {'hip_flexion_r','hip_adduction_r','knee_angle_r',...
    'ankle_angle_r','subtalar_angle_r','mtp_angle_r'};
vartype = 'Qs';

[fig_hand_1] = plot_figure_generic(R,model_info,coord_muscle_names_sel,vartype);


