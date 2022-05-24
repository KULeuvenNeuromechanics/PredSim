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
%       > dofs: cell array with coordinate names or muscle names
%       > variables: what to plot e.g. 'Qs', 'T_ID', 'lMtilde'
%       > savepath: path + name to save figure, no file extension
%       > filetyp: file extension: 'png', 'jpg', 'eps'
%
% 
% Original author: Lars D'Hondt
% Original date: 20/May/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

% generate colours, in this case a rainbow
colors = hsv(length(result_paths));

% pre-allocate figures to simplify code below
for j=1:length(figure_settings)
    fig_hands{j} = figure;
    fig_hands{j}.Name = figure_settings(j).name;
end

% loop over results
for i=1:length(result_paths)
    % load selected result
    load(result_paths{i},'R','model_info');

    % loop over figures
    for j=1:length(figure_settings)
        % check figure type
        if strcmp(figure_settings(j).dofs,'custom')
            % use custom figure function
            if strcmp(figure_settings(j).variables,'GRF')
                % use "plot_figure_grf.m"
                fig_hands{j} = plot_figure_grf(R,legend_names{i},colors(i,:),fig_hands{j});

            elseif strcmp(figure_settings(j).variables,'my_first_figure')
                % call you custom figure function here

            else
                % show template
                if i==1
                    fig_hands{j} = plot_figure_template(R,legend_names{i},colors(i,:),fig_hands{j});
                end
            end

        else % use generic figure function
            % replace "all_coords" by all coordinate names
            if strcmp(figure_settings(j).dofs,'all_coords')
                figure_settings(j).dofs = R.colheaders.coordinates;
            end
            % call function to plot generic figures
            fig_hands{j} = plot_figure_generic(R,model_info,figure_settings(j).dofs,...
                figure_settings(j).variables,legend_names{i},colors(i,:),fig_hands{j});
        end
    end


end


%% save figures if wanted
for j=1:length(figure_settings)
    if ~isempty(figure_settings(j).savepath) && ~isempty(figure_settings(j).filetype)
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



