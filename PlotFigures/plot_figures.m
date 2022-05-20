function [] = plot_figures(result_paths,legend_names,makeplot)
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
%   - makeplot -
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

colors = hsv(length(result_paths));

for i=1:length(result_paths)

    load(result_paths{i},'R','model_info');

    for j=1:length(makeplot)
        if strcmp(makeplot(j).dofs,'all_coords')
            makeplot(j).dofs = R.colheaders.coordinates;
        end
        
        if i==1
            fig_hands{j} = plot_figure_generic(R,model_info,makeplot(j).dofs,makeplot(j).variables,...
                legend_names{i},colors(i,:));

            fig_hands{j}.Name = makeplot(j).name;

        else

            fig_hands{j} = plot_figure_generic(R,model_info,makeplot(j).dofs,makeplot(j).variables,...
                legend_names{i},colors(i,:),fig_hands{j});

        end

    end


end


%% save figures if wanted
for j=1:length(makeplot)
    if ~isempty(makeplot(j).savepath) && ~isempty(makeplot(j).filetype)
        set(fig_hands{j},'PaperPositionMode','auto')
        for i=1:length(makeplot(j).filetype)
            if ~iscell(makeplot(j).filetype)
                filetype = {makeplot(j).filetype};
            else
                filetype = makeplot(j).filetype;
            end
            if strcmp(filetype(i),'png')
                print(fig_hands{j},makeplot(j).savepath,'-dpng','-r0')
            elseif strcmp(filetype(i),'jpg') || strcmp(filetype(i),'jpeg')
                print(fig_hands{j},makeplot(j).savepath,'-djpeg','-r0')
            elseif strcmp(filetype(i),'eps')
                print(fig_hands{j},makeplot(j).savepath,'-depsc')
            end
        end
    end
end