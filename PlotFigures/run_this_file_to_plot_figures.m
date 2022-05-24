
clear
close all
clc

[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
[pathRepoFolder,~,~] = fileparts(pathRepo);

%% General settings
% These settings will apply to all figures

% Construct a cell array with full paths to files with saved results for
% which you want to appear on the plotted figures.
results_folder = fullfile(pathRepoFolder,'PredSimResults');
result_paths{1} = fullfile([results_folder '\Fal_s1_mtp\Fal_s1_mtp_job148.mat']);
result_paths{2} = fullfile([results_folder '\Fal_s1_mtp_v2\Fal_s1_mtp_v2_job149.mat']);

% Cell array with legend name for each result
legend_names = {'test with old .dll', 'test with generated .dll'};

% Path to the folder where figures are saved
figure_folder = results_folder;

% Common part of the filename for all saved figures
figure_savename = 'test_dll';

%% Settings for each figure to be made
% "figure_settings" is a cell array where each cell contains a struct with
% the settings for a single figure.
% These settings are defined by several fields:
%   - name -
%   * String. Name assigned to the figure, and by default
%       appended to the filename when saving the figure.
%
%   - dofs -
%   * Cell array of strings. 
%

% 'C:\Users\u0150099\Documents\PredSim_dev\test_results\test_dll_angles'

% initilise the counter for dynamic indexing
fig_count = 1;

figure_settings(fig_count).name = 'all_angles';
figure_settings(fig_count).dofs = {'all_coords'};
figure_settings(fig_count).variables = {'Qs'};
figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
figure_settings(fig_count).filetype = {};
fig_count = fig_count+1;

% figure_settings(fig_count).name = 'selected_angles';
% figure_settings(fig_count).dofs = {'hip_flexion_r','hip_adduction_r','knee_angle_r',...
%     'ankle_angle_r','subtalar_angle_r','mtp_angle_r'};
% figure_settings(fig_count).variables = {'Qs'};
% figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
% figure_settings(fig_count).filetype = {};
% fig_count = fig_count+1;

% figure_settings(fig_count).name = 'torques';
% figure_settings(fig_count).dofs = {'all_coords'};
% figure_settings(fig_count).variables = {'T_ID'};
% figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
% figure_settings(fig_count).filetype = {};
% fig_count = fig_count+1;

figure_settings(fig_count).name = 'ankle_muscles';
figure_settings(fig_count).dofs = {'soleus_r','med_gas_r','lat_gas_r','tib_ant_r'};
figure_settings(fig_count).variables = {'a','FT','lMtilde','Wdot','Edot_gait'};
figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
figure_settings(fig_count).filetype = {};
fig_count = fig_count+1;

% figure_settings(fig_count).name = 'grfs';
% figure_settings(fig_count).dofs = {'custom'};
% figure_settings(fig_count).variables = {'GRF'};
% figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
% figure_settings(fig_count).filetype = {};
% fig_count = fig_count+1;

% figure_settings(fig_count).name = 'template';
% figure_settings(fig_count).dofs = {'custom'};
% figure_settings(fig_count).variables = {' '};
% figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
% figure_settings(fig_count).filetype = {};
% fig_count = fig_count+1;

%%

plot_figures(result_paths,legend_names,figure_settings);
