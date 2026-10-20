% test_plot_left_gait_cycle
%   Test script for plotting the left gait cycle with plot_figures.
%   Run the script cell by cell (Ctrl+Enter), first the settings cell,
%   then one of the options. Each option says what to check.
%
%   Option 1: right and left gait cycle in separate figures
%             (figure_settings.side = 'right' / 'left')
%   Option 2: default, no field "side" -> right gait cycle only
%   Option 3: left gait cycle added explicitly with add_left_gait_cycle
%             (R.left precomputed), right and left in the same figure
%             (plot functions called directly, not via plot_figures)
%
%   In options 1 and 2 the result has no R.left, so plot_figures computes
%   the left gait cycle.
%
% Original author: Ellis Van Can
% Original date: October 20, 2026

%% Settings (run this cell first)
clear
close all
clc

% Get paths for later use
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
[pathRepoFolder,~,~] = fileparts(pathRepo);

% Results to plot
results_folder = fullfile(pathRepoFolder,'PredSimResults');
result_paths{1} = fullfile(pathRepo,'Tests','ReferenceResults','Falisse_et_al_2022','Falisse_et_al_2022_paper.mat');
% result_paths{2} = fullfile(results_folder,'Falisse_et_al_2022_Results','Falisse_et_al_2022_v1.mat');

% Legend name for each result (no underscores: the legend interprets tex)
legend_names = {'Reference'};

% Saving figures
figure_folder = results_folder;
figure_savename = 'LeftGaitCycleTest';
filetype = {};   % e.g. {'png'} to save the figures

% Figures to make (same figures as run_this_file_to_plot_figures)
% dofs: full coordinate or muscle names (with _r and _l), or
%   'all_coords' -> all coordinates
%   'muscles_r'  -> all muscles of the right leg
%   'muscles_l'  -> all muscles of the left leg
%   'custom'     -> custom figure (here: 'GRF')
fig_count = 1;

base_settings(fig_count).name = 'all_angles';
base_settings(fig_count).dofs = {'all_coords'};
base_settings(fig_count).variables = {'Qs'};
fig_count = fig_count+1;

base_settings(fig_count).name = 'all_velocities';
base_settings(fig_count).dofs = {'all_coords'};
base_settings(fig_count).variables = {'Qdots'};
fig_count = fig_count+1;

base_settings(fig_count).name = 'all_accelerations';
base_settings(fig_count).dofs = {'all_coords'};
base_settings(fig_count).variables = {'Qddots'};
fig_count = fig_count+1;

base_settings(fig_count).name = 'all_activations_r';
base_settings(fig_count).dofs = {'muscles_r'};
base_settings(fig_count).variables = {'a'};
fig_count = fig_count+1;

base_settings(fig_count).name = 'all_activations_l';
base_settings(fig_count).dofs = {'muscles_l'};
base_settings(fig_count).variables = {'a'};
fig_count = fig_count+1;

base_settings(fig_count).name = 'selected_angles';
base_settings(fig_count).dofs = {'hip_flexion_r','hip_flexion_l','hip_adduction_r','hip_adduction_l',...
    'hip_rotation_r','hip_rotation_l','knee_angle_r','knee_angle_l',...
    'ankle_angle_r','ankle_angle_l','subtalar_angle_r','subtalar_angle_l',...
    'mtp_angle_r','mtp_angle_l'};
base_settings(fig_count).variables = {'Qs'};
fig_count = fig_count+1;

base_settings(fig_count).name = 'torques';
base_settings(fig_count).dofs = {'all_coords'};
base_settings(fig_count).variables = {'T_ID'};
fig_count = fig_count+1;

base_settings(fig_count).name = 'ankle_muscles';
base_settings(fig_count).dofs = {'soleus_r','soleus_l','med_gas_r','med_gas_l',...
    'lat_gas_r','lat_gas_l','tib_ant_r','tib_ant_l'};
base_settings(fig_count).variables = {'a','FT','lMtilde','Wdot','Edot_gait'};
fig_count = fig_count+1;

base_settings(fig_count).name = 'grfs';
base_settings(fig_count).dofs = {'custom'};
base_settings(fig_count).variables = {'GRF'};
fig_count = fig_count+1;

% save settings, the same for all figures
[base_settings.filetype] = deal(filetype);

%% Option 1: right and left gait cycle in separate figures
% Check:
%   - every figure exists twice: <name>_right and <name>_left
%   - in the _left figures the curves are shifted by about half a gait
%     cycle compared to the _right figures
%   - for a symmetric gait: knee_angle_r in <name>_right looks like
%     knee_angle_l in <name>_left (same for GRFs: right <-> left)
%   - grfs_left: the centre of pressure is the mirror image of grfs_right
close all
clear figure_settings
sides = {'right','left'};
fig_count = 1;
for s = 1:numel(sides)
    for j = 1:length(base_settings)
        % add the fields to a copy first, so every element has the same
        % fields before it is placed in the struct array
        fig_j = base_settings(j);
        fig_j.name = [base_settings(j).name '_' sides{s}];
        fig_j.side = sides{s};
        fig_j.savepath = fullfile(figure_folder,[figure_savename '_' fig_j.name]);
        figure_settings(fig_count) = fig_j;
        fig_count = fig_count+1;
    end
end
plot_figures(result_paths,legend_names,figure_settings);

%% Option 2: default, no field "side"
% Check:
%   - plot_figures runs without the field "side"
%   - the figures are the right gait cycle, identical to the <name>_right
%     figures of option 1 and to run_this_file_to_plot_figures
close all
clear figure_settings
figure_settings = base_settings;
for j = 1:length(figure_settings)
    figure_settings(j).savepath = fullfile(figure_folder,...
        [figure_savename '_' figure_settings(j).name '_default']);
end
assert(~isfield(figure_settings,'side'), 'figure_settings should not have the field side here.')
plot_figures(result_paths,legend_names,figure_settings);

%% Option 3: add left gait cycle explicitly, right and left in one figure
% The left gait cycle is added to R with add_left_gait_cycle. Each figure
% is plotted twice on the same figure handle: first with R (right gait
% cycle), then with R_left (left gait cycle). plot_figures is not used
% here; the plot functions are called directly.
% Check:
%   - R.left exists after add_left_gait_cycle
%   - every subplot shows two curves: right gait cycle and left gait cycle
%   - both curves have the same shape but are shifted by about half a
%     gait cycle
%   - grfs: the centre of pressure of the left cycle is the mirror image
%     of the right cycle, without missing parts
close all

% create figures
fig_hands = cell(1,length(base_settings));
for j = 1:length(base_settings)
    fig_hands{j} = figure;
    fig_hands{j}.Name = [base_settings(j).name '_both'];
end

% two colours per result: one for the right cycle, one for the left cycle
colors = lines(2*length(result_paths));

for i = 1:length(result_paths)
    load(result_paths{i},'R','model_info');

    % add left gait cycle explicitly
    R = add_left_gait_cycle(R,model_info);
    assert(isfield(R,'left'), 'add_left_gait_cycle did not add R.left.')

    % R_left: copy of R with the fields of the left gait cycle (R.left)
    R_left = rmfield(R,'left');
    categories = fieldnames(R.left);
    for c = 1:numel(categories)
        cat_c = categories{c};
        if ~isstruct(R.left.(cat_c))
            continue % e.g. idx_GC
        end
        fn = fieldnames(R.left.(cat_c));
        for f = 1:numel(fn)
            val = R.left.(cat_c).(fn{f});
            if isstruct(val)
                % e.g. R.left.metabolics.Bhargava2004
                fn_sub = fieldnames(val);
                for k = 1:numel(fn_sub)
                    R_left.(cat_c).(fn{f}).(fn_sub{k}) = val.(fn_sub{k});
                end
            else
                R_left.(cat_c).(fn{f}) = val;
            end
        end
    end

    leg_right = [legend_names{i} ' - right gait cycle'];
    leg_left  = [legend_names{i} ' - left gait cycle'];

    for j = 1:length(base_settings)
        dofs = base_settings(j).dofs;

        % GRF figure
        if any(strcmp(dofs,'custom'))
            if any(strcmp(base_settings(j).variables,'GRF'))
                fig_hands{j} = plot_figure_grf(R,leg_right,colors(2*i-1,:),fig_hands{j});
                fig_hands{j} = plot_figure_grf(R_left,leg_left,colors(2*i,:),fig_hands{j});
            end
            continue
        end

        % replace keywords by the full list of names
        if any(strcmp(dofs,'all_coords'))
            dofs = R.colheaders.coordinates;
        elseif any(strcmp(dofs,'muscles_r')) || any(strcmp(dofs,'muscles_l'))
            if any(strcmp(dofs,'muscles_r'))
                suffix = '_r';
            else
                suffix = '_l';
            end
            muscle_names = R.colheaders.muscles;
            is_side = false(length(muscle_names),1);
            for im = 1:length(muscle_names)
                is_side(im) = strcmp(muscle_names{im}(end-1:end),suffix);
            end
            dofs = muscle_names(is_side);
        end

        % right gait cycle, then left gait cycle in the same subplots
        fig_hands{j} = plot_figure_generic(R,model_info,dofs,...
            base_settings(j).variables,leg_right,colors(2*i-1,:),fig_hands{j});
        fig_hands{j} = plot_figure_generic(R_left,model_info,dofs,...
            base_settings(j).variables,leg_left,colors(2*i,:),fig_hands{j});
    end
end

% save figures if wanted
for j = 1:length(base_settings)
    for k = 1:length(filetype)
        exportgraphics(fig_hands{j},fullfile(figure_folder,...
            [figure_savename '_' base_settings(j).name '_both.' filetype{k}]),'Resolution',300);
    end
end