
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
result_paths{1} = fullfile([results_folder '\CP3\CP3_v4.mat']);
result_paths{2} = fullfile([results_folder '\CP3\CP3_v5.mat']);
result_paths{3} = fullfile([results_folder '\CP3\CP3_v6.mat']);
result_paths{4} = fullfile([results_folder '\CP3\CP3_v7.mat']);

% Cell array with legend name for each result
legend_names = {'half gait cycle (Nmesh = 50)','full gait cycle (Nmesh = 50)',...
    'full gait cycle (Nmesh = 100)','full gait cycle (Nmesh = 100, QR)'};

% Path to the folder where figures are saved
figure_folder = results_folder;

% Common part of the filename for all saved figures
figure_savename = 'CP3_MRI2_comparemeshes_IG';

%% Settings for each figure to be made
% "figure_settings" is a cell array where each cell contains a struct with
% the settings for a single figure.
% These settings are defined by several fields:
%   - name -
%   * String. Name assigned to the figure, and by default
%   appended to the filename when saving the figure.
%
%   - dofs -
%   * Cell array of strings. Can contain coordinate names OR muscle names.
%   Alternatively, 'all_coords' will use all coordinates from the 1st
%   result. Enter 'custom' to use variables that do not exist for individual
%   coordinates or muscles.
%
%   - variables -
%   * Cell array of strings. Containsone or more variable names. e.g. 'Qs'
%   to plot coordinate positions, 'a' to plot muscle activity. Variables
%   that do not rely on coordinates or muscles (e.g. GRFs)
%
%   - savepath -
%   * String. Full path + filename used to save the figure. Does not
%   include file extension.
%
%   - filetype -
%   * Cell array of strings. File extensions to save the figure as, leave
%   empty to not save the figure. Supported types are: 'png', 'jpg', 'eps'
%
%



% initilise the counter for dynamic indexing
fig_count = 1;

figure_settings(fig_count).name = 'all_angles';
figure_settings(fig_count).dofs = {'all_coords'};
figure_settings(fig_count).variables = {'Qs'};
figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
figure_settings(fig_count).filetype = {'jpeg'};
fig_count = fig_count+1;

% figure_settings(fig_count).name = 'selected_angles';
% figure_settings(fig_count).dofs = {'hip_flexion_r','hip_adduction_r','hip_rotation_r','knee_angle_r',...
%     'ankle_angle_r','subtalar_angle_r','mtp_angle_r'};
% figure_settings(fig_count).variables = {'Qs'};
% figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
% figure_settings(fig_count).filetype = {};
% fig_count = fig_count+1;

figure_settings(fig_count).name = 'torques';
figure_settings(fig_count).dofs = {'all_coords'};
figure_settings(fig_count).variables = {'T_ID'};
figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
figure_settings(fig_count).filetype = {'jpeg'};
fig_count = fig_count+1;

figure_settings(fig_count).name = 'ankle_muscles';
figure_settings(fig_count).dofs = {'soleus_r','med_gas_r','lat_gas_r','tib_ant_r'};
figure_settings(fig_count).variables = {'a','FT','lMtilde','Wdot','Edot_gait'};
figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
figure_settings(fig_count).filetype = {'jpeg'};
fig_count = fig_count+1;

figure_settings(fig_count).name = 'evert_muscles_l';
figure_settings(fig_count).dofs = {'ext_dig_l','per_brev_l','per_long_l','per_tert_l'};
figure_settings(fig_count).variables = {'a','FT','lMtilde','Wdot','Edot_gait'};
figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
figure_settings(fig_count).filetype = {'jpeg'};
fig_count = fig_count+1;

figure_settings(fig_count).name = 'grfs';
figure_settings(fig_count).dofs = {'custom'};
figure_settings(fig_count).variables = {'GRF'};
figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
figure_settings(fig_count).filetype = {'jpeg'};
fig_count = fig_count+1;

% figure_settings(fig_count).name = 'template';
% figure_settings(fig_count).dofs = {'custom'};
% figure_settings(fig_count).variables = {' '};
% figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
% figure_settings(fig_count).filetype = {};
% fig_count = fig_count+1;

%%

plot_figures(result_paths,legend_names,figure_settings);
