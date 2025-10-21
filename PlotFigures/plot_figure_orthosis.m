function [fig_hand] = plot_figure_orthosis(R, namesToPlot, varargin)
% --------------------------------------------------------------------------
% plot_figure_template
%   Plots orthosis states, controls, and/or exoskeleton joint torques
% 
% INPUT:
%   - R -
%   * struct with simulation results
%
%   - namesToPlot -
%   * cell array of strings (names of orthosis states, controls, or DOFs)
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
% Original date: 24/May/2022
%
% Last edit by: ChatGPT (GPT-5)
% Last edit date: 21/Oct/2025
% --------------------------------------------------------------------------

%% Early exit if no orthosis data
if ~isfield(R, 'orthosis') || isempty(R.orthosis)
    warning('R.orthosis does not exist. Exiting plot_figure_orthosis without plotting.');
    fig_hand = []; % return empty handle
    return
end

%% Default settings
colr = [];
legName = R.S.misc.result_filename;
legInt = 'none';

%% Input parsing
for i = 1:length(varargin)
    if isa(varargin{i},'matlab.ui.Figure')
        fig_hand = varargin{i};
    elseif isa(varargin{i},'double') && length(varargin{i})==3
        colr = varargin{i};
    elseif isa(varargin{i},'string') || isa(varargin{i},'char')
        legName = varargin{i};
        legInt = 'tex';
    end
end

%% Determine figure
if ~exist('fig_hand','var') || isempty(fig_hand)
    fig_hand = figure('Name','Orthosis signals');
end

% Create axes explicitly
ax = axes(fig_hand); 
hold(ax,'on');


% Use number of mesh points from solver instead of relying on orthosis.states (may not exist)
percentage_gait_cycle = linspace(0,100,R.S.solver.N_meshes);

% available names
if isfield(R, 'S') && isfield(R.S, 'orthosis')
    if isfield(R.S.orthosis, 'stateNames_all')
        stateNames = R.S.orthosis.stateNames_all;
    else
        stateNames = {};
    end
    if isfield(R.S.orthosis, 'controlNames_all')
        controlNames = R.S.orthosis.controlNames_all;
    else
        controlNames = {};
    end
else
    stateNames = {};
    controlNames = {};
end

if isfield(R, 'colheaders') && isfield(R.colheaders, 'coordinates')
    coordNames = R.colheaders.coordinates; % used for T_coord matching
else
    coordNames = {};
end

%% Prepare default color cycling
defaultColors = get(gca,'ColorOrder');
nColors = size(defaultColors,1);

%% Loop through requested names
for i = 1:length(namesToPlot)
    name = namesToPlot{i};
    data = [];
    src  = '';

    % --- match in orthosis states
    idx_state = find(strcmp(name, stateNames));
    if ~isempty(idx_state)
        data = R.orthosis.states(:, idx_state);
        src = 'state';
    end

    % --- match in orthosis controls
    if isempty(data)
        idx_ctrl = find(strcmp(name, controlNames));
        if ~isempty(idx_ctrl)
            data = R.orthosis.controls(:, idx_ctrl);
            src = 'control';
        end
    end

    % --- match in combined torques (T_coord)
    if isempty(data) && isfield(R.orthosis,'combined') && ...
            isfield(R.orthosis.combined,'T_coord')
        % like plot_figure_generic, look for exact match in coordinate headers
        idx_torque = find(strcmp(R.colheaders.coordinates, name));
        if ~isempty(idx_torque)
            data = R.orthosis.combined.T_coord(:, idx_torque);
            src = 'exo torque';
        end
    end

    % --- handle not found
    if isempty(data)
        warning('Name "%s" not found in states, controls, or DOFs.', name);
        continue;
    end

    % Determine color
    if isempty(colr)
        colorToUse = defaultColors(mod(i-1,nColors)+1, :);
    else
        colorToUse = colr;
    end

    % --- plot data
    plot(percentage_gait_cycle, data, ...
        'DisplayName', sprintf('%s (%s)', name, src), ...
        'Color', colorToUse, 'LineWidth', 1.5);
end

%% Format axes
xlabel('Gait cycle (%)')
ylabel('Value')
title(legName, 'Interpreter','none')
legend('Location','best','Interpreter',legInt)
grid on
box on

end








