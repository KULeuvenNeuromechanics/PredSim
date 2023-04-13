
%% Clear workspace
if exist('result_paths','var')
    clearvars -Except result_paths
else
    clear
end
close all
clc

%% Get paths for later use
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
[pathRepoFolder,~,~] = fileparts(pathRepo);

%% General settings
% These settings will apply to all figures

% Construct a cell array with full paths to files with saved results for
% which you want to appear on the plotted figures.
results_folder = fullfile(pathRepoFolder,'PredSimResults');
% result_paths{1} = fullfile(pathRepo,'Tests','Falisse_et_al_2022_Results','Falisse_et_al_2022_v1.mat');
result_paths{1} = fullfile(results_folder,'DHondt_2023_3seg\DHondt_2023_3seg_ref.mat');
% result_paths{end+1} = fullfile(results_folder,'DHondt_2023_3seg\DHondt_2023_3seg_job649.mat');
% result_paths{end+1} = fullfile(results_folder,'DHondt_2023_3seg\DHondt_2023_3seg_job650.mat');
% result_paths{end+1} = fullfile(results_folder,'DHondt_2023_3seg\DHondt_2023_3seg_job651.mat');
% result_paths{end+1} = fullfile(results_folder,'DHondt_2023_3seg\DHondt_2023_3seg_job653.mat');
% result_paths{end+1} = fullfile(results_folder,'DHondt_2023_3seg\DHondt_2023_3seg_job654.mat');
% result_paths{end+1} = fullfile(results_folder,'DHondt_2023_3seg\DHondt_2023_3seg_job661.mat');
% result_paths{end+1} = fullfile(results_folder,'DHondt_2023_3seg\DHondt_2023_3seg_job663.mat');
% result_paths{end+1} = fullfile(results_folder,'DHondt_2023_3seg\DHondt_2023_3seg_job665.mat');
result_paths{end+1} = fullfile(results_folder,'DHondt_2023_3seg\DHondt_2023_3seg_job667.mat');


% Cell array with legend name for each result
legend_names = {'reference'};

% Path to the folder where figures are saved
figure_folder = results_folder;

% Common part of the filename for all saved figures
figure_savename = 'MyFirstPredictiveSimulationFigure';

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
figure_settings(fig_count).filetype = {};
fig_count = fig_count+1;

% figure_settings(fig_count).name = 'all_angles';
% figure_settings(fig_count).dofs = {'all_coords'};
% figure_settings(fig_count).variables = {'Qdots'};
% figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
% figure_settings(fig_count).filetype = {};
% fig_count = fig_count+1;

% figure_settings(fig_count).name = 'all_angles';
% figure_settings(fig_count).dofs = {'all_coords'};
% figure_settings(fig_count).variables = {'Qddots'};
% figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
% figure_settings(fig_count).filetype = {};
% fig_count = fig_count+1;

figure_settings(fig_count).name = 'all_activations';
figure_settings(fig_count).dofs = {'muscles_r'};
figure_settings(fig_count).variables = {'a'};
figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
figure_settings(fig_count).filetype = {};
fig_count = fig_count+1;

% figure_settings(fig_count).name = 'selected_angles';
% figure_settings(fig_count).dofs = {'hip_flexion_r','hip_adduction_r','hip_rotation_r','knee_angle_r',...
%     'ankle_angle_r','subtalar_angle_r','mtp_angle_r'};
% figure_settings(fig_count).variables = {'Qs'};
% figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
% figure_settings(fig_count).filetype = {'jpeg'};
% fig_count = fig_count+1;

figure_settings(fig_count).name = 'torques';
figure_settings(fig_count).dofs = {'all_coords'};
figure_settings(fig_count).variables = {'T_ID'};
figure_settings(fig_count).savepath = fullfile(figure_folder,[figure_savename '_' figure_settings(fig_count).name]);
figure_settings(fig_count).filetype = {};
fig_count = fig_count+1;

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
