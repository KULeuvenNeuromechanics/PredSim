function [fig_hand] = plot_figure_generic(R,model_info,coord_muscle_names_sel,vartype,varargin)
% --------------------------------------------------------------------------
% plot_figure_generic
%   Plots a figure with any combination of coordinate-related variables, or
%   muscle-related variables. This function is very generic, you are
%   adviced to create your own figure template function rather than
%   customising this one.
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
%       coordinates take priority in case of mixed input
%
%   - vartype - 
%   * cell array with names of variables, or name of a single variable to 
%   plot for selected coordinates (e.g. Qs to plot angles)
%
%   - optional inputs - 
%   * These inputs are recognised by variable type, so their order does not matter
%   > figure handle: add plots to an existing figure (default = new figure)
%   > array of 3 doubles: use as RGB triplet for the new plot (default = blue)
%   > string or character array: use as name in the legend (default = filename)
%   > 'right': align figure window to right side of screen (default = left)
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
% Default settings:
% figure on left side of screen 
fig_pos = 0;
% default color is blue
colr = [0 0.4470 0.7410];
% legend name
legName = R.S.post_process.result_filename;
% use no interpreter for legend
lgInt = 'none';

% Use optional inputs to overwrite default settings
for i=1:length(varargin)
    if isa(varargin{i},'matlab.ui.Figure')
        % figure handle
        fig_hand = varargin{i};

    elseif isa(varargin{i},'double') && length(varargin{i})==3
        % RGB triplet
        colr = varargin{i};

    elseif isa(varargin{i},'char') && strcmp(varargin{i},'right')
        % figure on right side of screen
        fig_pos = 1;

    elseif isa(varargin{i},'string') || isa(varargin{i},'char')
        % use given name for legend
        legName = varargin{i};
        % interpret tex
        lgInt = 'tex';

    end
end

% make sure vartype is a cell array
if isa(vartype,'string') || isa(vartype,'char')
    vartype1 = {vartype};
    vartype = vartype1;
end

%% find results
nvar = length(vartype);
% preallocate cell array to place data for plotting
ydata = cell(nvar,1);

% when plotting coordinates...
if ~isempty(intersect(R.colheaders.coordinates,coord_muscle_names_sel))
    colheaders = R.colheaders.coordinates;
    for i=1:nvar
        %... expect kinematics and kinetics
        if isfield(R.kinematics,vartype{i})
            ydata{i} = R.kinematics.(vartype{i});
        elseif isfield(R.kinetics,vartype{i})
            ydata{i} = R.kinetics.(vartype{i});
        end
    end

% when plotting muscles...
elseif ~isempty(intersect(R.colheaders.muscles,coord_muscle_names_sel))
    colheaders = R.colheaders.muscles;
    %... select a metabolic energy model used for simulation
    if isfield(R.metabolics,R.S.metabolicE.model)
        metabolic_model = R.S.metabolicE.model;
    else
        metabolic_model = 'Bhargava2004';
    end
    for i=1:nvar
        %... expect muscles and metabolics
        if isfield(R.muscles,vartype{i})
            ydata{i} = R.muscles.(vartype{i});
        elseif isfield(R.metabolics.(metabolic_model),vartype{i})
            ydata{i} = R.metabolics.(metabolic_model).(vartype{i});
        end
    end
end

% Remove vartype elements that are not recognised as a variable
% idx_found = find(~cellfun('isempty',ydata));
% vartype = vartype(idx_found);
% ydata = ydata(idx_found);

% Remove coordinate or muscle names that were not recognised
% coord_muscle_names_sel = intersect(colheaders,coord_muscle_names_sel,'stable');

%% Determine figure layout
% number of horizontal subplots: column for each coordinate/muscle
nh = length(coord_muscle_names_sel);
% number of vertical subplots: row for every vartype
nv = length(vartype);
% If there would be only 1 row, let this variable run through multiple rows
if nv == 1
    nh1 = ceil(sqrt(nh));
    nv = ceil(nh/nh1);
    nh = nh1;
end

%% Set figure size and position on screen
% fraction of screen size
fsh = min(nh/8,1);
fsv = min(nv/5,1);
% screen size
scs = get(0,'ScreenSize');
scs(4) = scs(4) - 150;
scs(3) = scs(3) - 4;
% figure size
fig_size = [fsh*scs(3), fsv*scs(4)];
% figure position
fposh = fig_pos*(scs(3)-fig_size(1))+2;
fposv = (scs(4)-fig_size(2))/2;
fig_pos = [fposh,fposv];

%% Create or select figure
% If there was no figure handle provided, create a figure for this handle
if ~exist('fig_hand','var')
    fig_hand = figure;
    bool_1st = 1;
else
    bool_1st = 0;
end
% set position and size of figure
fig_hand.Position = [fig_pos,fig_size];
% select figure
figure(fig_hand)

%%
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

        if cnt_y == 1
            title(coord_muscle_names_sel(cnt_var),'Interpreter','none')
        end

    end

%     if cnt_y == 1
%         title(coord_muscle_names_sel(cnt_var),'Interpreter','none')
%     end

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

lgh = legend('Location','northwest','Interpreter',lgInt);

if nsp<nv*nh
    lgh.Position(1) = lgh.Position(1)+1/(nh+1);
else
%     lgh.Position(2) = lgh.Position(2)-1/(nv*2);
end

