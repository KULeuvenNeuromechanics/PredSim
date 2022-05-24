function [fig_hand] = plot_figure_template(R,varargin)
% --------------------------------------------------------------------------
% plot_figure_template
%   Please use this
% 
% INPUT:
%   - R -
%   * struct with simulation results
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
% Last edit by:
% Last edit date:
% --------------------------------------------------------------------------


%% Default settings
% default color is blue
colr = [0 0.4470 0.7410];
% default legend name is filename
legName = R.S.post_process.result_filename;
% by default, use no interpreter for legend, because default name contains "_"
legInt = 'none';

%% Input settings
% loop through optional inputs
for i=1:length(varargin)
    if isa(varargin{i},'matlab.ui.Figure')
        % use existing figure if handle is provided
        fig_hand = varargin{i};
    elseif isa(varargin{i},'double') && length(varargin{i})==3
        % use provided
        colr = varargin{i};
    elseif isa(varargin{i},'string') || isa(varargin{i},'char')
        legName = varargin{i};
        legInt = 'tex';
    end
end

%% Set figure size and position on screen
% screen size of your main screen, in pixels (3rd element is width, 4th element is height)
scs = get(0,'ScreenSize');
% reduce height to create margin for header
scs(4) = scs(4) - 150;
% height and width of figure are half of respective screen dimensions
fig_size = [scs(3)/2,scs(4)/2];
% bottom left corner of figure is as quarter of screen dimensions
fig_pos = [scs(3)/4,scs(4)/4];


%% Create or select figure
% If there was no figure handle provided, create a figure for this handle
if ~exist('fig_hand','var')
    fig_hand = figure;
end
% set position and size of figure
fig_hand.Position = [fig_pos,fig_size];
% select figure
figure(fig_hand)


%% Plotting
% example plot
percentage_gait_cycle = linspace(1,100,size(R.ground_reaction.GRF_r,1));

subplot(1,2,1)
hold on
dist_forward = R.kinematics.Qs(:,strcmp(R.colheaders.coordinates,'pelvis_tx'));
plot(percentage_gait_cycle,dist_forward,'color',colr,'DisplayName',legName)
title('Forward distance of floating base','Interpreter','none')
ylabel('Distance (m)')
xlabel('Gait cycle (%)')
legend('Location','best','Interpreter',legInt)


subplot(1,2,2)
text(0.05,0.5,...
    {'You can find a template to create additional figure','set-ups in /PlotFigures/plot_figure_template.m'},...
    'Interpreter','none');







