function [] = plot_figures(result_paths,legend_names,figure_settings)
% --------------------------------------------------------------------------
% plot_figures
%   This functions plots the figures specified in makeplot, for all results
%   given in result_paths.
% 
% INPUT:
%   - result_paths -
%   * cell array with paths to result files
%
%   - legend_names -
%   * cell array with legend name for each result
% 
%   - figure_settings -
%   * cell array with structs that each describe a figure
%       > name: name of the figure
%       > dofs: cell array with specific coordinate names or muscle names
%           'all_coords', 'muscles_r' or 'muscles_l' are also accepted   
%       > variables: what to plot e.g. 'Qs', 'T_ID', 'lMtilde'
%       > side: (optional) 'right' (default) or 'left'. Selects which
%           gait cycle is plotted. 'left' uses R.left (see
%           add_left_gait_cycle.m); if a result has no R.left, it is
%           computed within this function
%       > savepath: path + name to save figure, no file extension
%       > filetyp: file extension: 'png', 'jpg', 'eps'
%
% 
% Original author: Lars D'Hondt
% Original date: 20/May/2022
%
% Last edit by: Ellis Van Can
% Last edit date: October 20, 2026
% --------------------------------------------------------------------------
% This file is part of PredSim.
% 
% PredSim: A Framework for Rapid Predictive Simulations of Locomotion
% Copyright (c) 2026 KU Leuven
% 
% PredSim is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Affero General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% 
% PredSim is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public 
% License for more details.
% 
% You should have received a copy of the GNU Affero General Public License 
% along with PredSim. If not, see <https://www.gnu.org/licenses/>.
% --------------------------------------------------------------------------

% generate colours, in this case a rainbow
colors = hsv(length(result_paths));

% default side is right
if ~isfield(figure_settings,'side')
    [figure_settings.side] = deal('right');
    
end

% pre-allocate figures to simplify code below
for j=1:length(figure_settings)
    fig_hands{j} = figure;
    fig_hands{j}.Name = figure_settings(j).name;
end

% loop over results
for i=1:length(result_paths)
    % load selected result
    load(result_paths{i},'R','model_info');

    if length(legend_names)<i
        legend_names{i} = replace(R.S.misc.result_filename,'_',' ');
    end

    % prepare left gait cycle once per result, only if a figure needs it
    % R_left is a copy of R in which the time-varying fields are replaced
    % by those of the left gait cycle (R.left), so all plot functions work
    % unchanged.
    R_left = [];
    if any(strcmpi({figure_settings.side},'left'))
        % older results (e.g. reference results) have no left cycle yet
        if ~isfield(R,'left')
            R = add_left_gait_cycle(R,model_info);
        end
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
                    % e.g. R.left.metabolics.Bhargava2004: overwrite per
                    % field, so non time-varying fields (e.g. COT) remain
                    fn_sub = fieldnames(val);
                    for k = 1:numel(fn_sub)
                        R_left.(cat_c).(fn{f}).(fn_sub{k}) = val.(fn_sub{k});
                    end
                else
                    R_left.(cat_c).(fn{f}) = val;
                end
            end
        end
    end
    
    % loop over figures
    for j=1:length(figure_settings)
        % select the gait cycle to plot
        if strcmpi(figure_settings(j).side,'left')
            R_plot = R_left;
        else
            R_plot = R;
        end

        % check figure type
        if strcmp(figure_settings(j).dofs,'custom')
            % use custom figure function
            if strcmp(figure_settings(j).variables,'GRF')
                % use "plot_figure_grf.m"
                fig_hands{j} = plot_figure_grf(R_plot,legend_names{i},colors(i,:),fig_hands{j});

            elseif strcmp(figure_settings(j).variables,'my_first_figure')
                % call you custom figure function here

            else
                % show template
                if i==1
                    fig_hands{j} = plot_figure_template(R_plot,legend_names{i},colors(i,:),fig_hands{j});
                end
            end

        else % use generic figure function
            % replace "all_coords" by all coordinate names
            if strcmp(figure_settings(j).dofs,'all_coords')
                figure_settings(j).dofs = R_plot.colheaders.coordinates;
            end
            % replace "muscles_r" / "muscles_l" by all muscle names of that side
            if any(strcmp(figure_settings(j).dofs,'muscles_r')) || ...
                    any(strcmp(figure_settings(j).dofs,'muscles_l'))
                if strcmp(figure_settings(j).dofs,'muscles_r')
                    suffix = '_r';
                else
                    suffix = '_l';
                end
                muscle_names = R_plot.colheaders.muscles;
                is_side = false(length(muscle_names),1);
                for im = 1:length(muscle_names)
                    is_side(im) = strcmp(muscle_names{im}(end-1:end),suffix);
                end
                figure_settings(j).dofs = muscle_names(is_side);
            end
            % call function to plot generic figures
            fig_hands{j} = plot_figure_generic(R_plot,model_info,figure_settings(j).dofs,...
                figure_settings(j).variables,legend_names{i},colors(i,:),fig_hands{j});
        end
    end

end


%% save figures if wanted
for j=1:length(figure_settings)
    if ~isempty(figure_settings(j).savepath) && ~isempty(figure_settings(j).filetype)
        % older versions of matlab do not have exportgraphics
        if exist('exportgraphics')
            for i=1:length(figure_settings(j).filetype)
                if ~iscell(figure_settings(j).filetype)
                    filetype = {figure_settings(j).filetype};
                else
                    filetype = figure_settings(j).filetype;
                end
                    exportgraphics(fig_hands{j},[figure_settings(j).savepath '.' filetype{i}],'Resolution',300);
            end

        else
            set(fig_hands{j},'PaperPositionMode','auto')
            for i=1:length(figure_settings(j).filetype)
                if ~iscell(figure_settings(j).filetype)
                    filetype = {figure_settings(j).filetype};
                else
                    filetype = figure_settings(j).filetype;
                end
                if strcmp(filetype(i),'png')
                    print(fig_hands{j},figure_settings(j).savepath,'-dpng','-r0')
                elseif strcmp(filetype(i),'jpg') || strcmp(filetype(i),'jpeg')
                    print(fig_hands{j},figure_settings(j).savepath,'-djpeg','-r0')
                elseif strcmp(filetype(i),'eps')
                    print(fig_hands{j},figure_settings(j).savepath,'-depsc')
                end
            end
        end
    end
end



