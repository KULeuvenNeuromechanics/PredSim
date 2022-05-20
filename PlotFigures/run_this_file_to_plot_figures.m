
clear
close all
clc

[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);

%% Select results to plot
results_folder = [pathRepo '\test_results'];
result_paths{1} = fullfile([results_folder '\Fal_s1_mtp\Fal_s1_mtp_job148.mat']);
result_paths{2} = fullfile([results_folder '\Fal_s1_mtp_v2\Fal_s1_mtp_v2_job149.mat']);

legend_names = {'test with old .dll', 'test with generated .dll'};

%% Save figures
figure_folder = results_folder;
figure_savename = 'test_dll';

%% Select plots to make
makeplot(1).name = 'angles';
makeplot(1).dofs = {'all_coords'};
makeplot(1).variables = {'Qs'};
makeplot(1).savepath = fullfile(figure_folder,[figure_savename '_' makeplot(1).name]);
makeplot(1).filetype = {};

makeplot(2).name = 'torques';
makeplot(2).dofs = {'all_coords'};
makeplot(2).variables = {'T_ID'};
makeplot(2).savepath = fullfile(figure_folder,[figure_savename '_' makeplot(2).name]);
makeplot(2).filetype = {};

% makeplot(3).name = 'ankle_muscles';
% makeplot(3).dofs = {'soleus_r','med_gas_r','lat_gas_r','tib_ant_r'};
% makeplot(3).variables = {'a','FT','lMtilde','Wdot','Edot_gait'};
% makeplot(3).savepath = fullfile(figure_folder,[figure_savename '_' makeplot(3).name]);
% makeplot(3).filetype = {};



%%

plot_figures(result_paths,legend_names,makeplot);
