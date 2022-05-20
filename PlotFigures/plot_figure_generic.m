function [fig_hand] = plot_figure_generic(R,model_info,coord_muscle_names_sel,vartype,varargin)
% --------------------------------------------------------------------------
% plot_figure_single_property
%   Plots a figure that shows a single property (e.g. 
% 
% INPUT:
%   - R -
%   * struct with simulation results
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
%   - coord_muscle_names_sel - 
%   * cell array with names of _coordinates OR muscles_ selected to plot
%
%   - vartype - 
%   * cell array with names of variables, or name of a single variable to 
%   plot for selected coordinates (e.g. Qs to plot angles)
%
%   - optional inputs - 
%   * pass a figure handle to add plots to an existing figure, otherwise
%       this function will make a new figure
%   * pass an array of 3 doubles to use as RGB triplet for the new plot
%   * pass a string to use as name in the legend
%
% OUTPUT:
%   - fig_hand -
%   * handle of the plotted figure
%
% Original author: Lars D'Hondt
% Original date: 20/May/2022
%
% Last edit by:
% Last edit date:
% --------------------------------------------------------------------------

%% settings
% figure on left side of screen 
fig_pos = 0;
% default color is blue
colr = [0 0.4470 0.7410];
% legend name
legName = R.S.post_process.result_filename;
% use no interpreter for legend
lgInt = 'none';

% identify and unpack varargin
for i=1:length(varargin)
    if isa(varargin{i},'matlab.ui.Figure')
        fig_hand = varargin{i};
    elseif isa(varargin{i},'double') && length(varargin{i})==3
        colr = varargin{i};
    elseif isa(varargin{i},'string') || isa(varargin{i},'char')
        legName = varargin{i};
        lgInt = 'tex';
    elseif strcmp(varargin{i},'right')
        fig_pos = 1;
    end
end

if isa(vartype,'string') || isa(vartype,'char')
    vartype1 = {vartype};
    vartype = vartype1;
end

% number of horizontal subplots
nh = length(coord_muscle_names_sel);
% number of vertical subplots
nv = length(vartype);

%% find results
ydata = cell(nv,1);
% when plotting coordinates
if ~isempty(intersect(R.colheaders.coordinates,coord_muscle_names_sel))
    colheaders = R.colheaders.coordinates;
    for i=1:nv
        if isfield(R.kinematics,vartype{i})
            ydata{i} = R.kinematics.(vartype{i});
        elseif isfield(R.kinetics,vartype{i})
            ydata{i} = R.kinetics.(vartype{i});
        end
    end

% when plotting muscles
elseif ~isempty(intersect(R.colheaders.muscles,coord_muscle_names_sel))
    colheaders = R.colheaders.muscles;
    for i=1:nv
        if isfield(R.muscles,vartype{i})
            ydata{i} = R.muscles.(vartype{i});
        elseif isfield(R.metabolics.Bhargava2004,vartype{i})
            ydata{i} = R.metabolics.Bhargava2004.(vartype{i});
        end
    end
end

idx_found = find(~cellfun('isempty',ydata));
vartype = vartype(idx_found);
ydata = ydata(idx_found);
nv = length(idx_found);

%% determine figure layout
% adapt subplot distribution
if nv == 1
    nh1 = ceil(sqrt(nh));
    nv = ceil(nh/nh1);
    nh = nh1;
end

% create figure
if ~exist('fig_hand','var')
    % fraction of screen size
    fsh = min(nh/8,1);
    fsv = min(nv/5,1);
    % screen size
    scs = get(0,'ScreenSize');
    scs(4) = scs(4) - 150;
    scs(3) = scs(3) - 4;
    % figure size
    fsize = [fsh*scs(3), fsv*scs(4)];
    % figure position
    fposh = fig_pos*(scs(3)-fsize(1))+2;
    fposv = (scs(4)-fsize(2))/2;
    % create
    fig_hand = figure('Position',[fposh,fposv,fsize]);
end


figure(fig_hand)
cnt_var = 1;
cnt_y = 1;

nsp = length(coord_muscle_names_sel)*length(vartype);

for i=1:nsp

    idx = find(strcmp(colheaders,coord_muscle_names_sel(cnt_var)));
    if ~isempty(idx)
        idx = idx(1);
        y_i = ydata{cnt_y}(:,idx);
        x_i = linspace(1,100,length(y_i));

        subplot(nv,nh,i)
        hold on
        plot(x_i,y_i,'color',colr,'DisplayName',legName)

    end

    if cnt_y == 1
        title(coord_muscle_names_sel(cnt_var),'Interpreter','none')
    end

    if mod(i,nh) == 1
        ylabel(vartype(cnt_y),'Interpreter','none')
    end

    if i > nsp-nh
        xlabel('Gait cycle (%)')
    end

    cnt_var = cnt_var + 1;
    if cnt_var > length(coord_muscle_names_sel)
        cnt_var = 1;
        cnt_y = cnt_y + 1;
    end

end

legend('Interpreter',lgInt);



