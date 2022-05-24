function [fig_hand] = plot_figure_grf(R,varargin)
% --------------------------------------------------------------------------
% plot_figure_grf
%   Plots a figure that shows total ground reaction forces and centre of
%   pressure.
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
% Original date: 23/May/2022
%
% Last edit by:
% Last edit date:
% --------------------------------------------------------------------------


%% settings
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

x = linspace(1,100,size(R.ground_reaction.GRF_r,1));

for i=1:3
    subplot(3,3,i)
    hold on
    plot(x,R.ground_reaction.GRF_r(:,i)/R.misc.body_weight*100,'color',colr,'DisplayName',legName)
    title([R.colheaders.GRF{i} ' right'],'Interpreter','none')
end

legend('Location','best','Interpreter',lgInt)

subplot(3,3,1)
ylabel('GRF (% BW)')

for i=1:3
    subplot(3,3,i+3)
    hold on
    plot(x,R.ground_reaction.GRF_l(:,i)/R.misc.body_weight*100,'color',colr,'DisplayName',legName)
    title([R.colheaders.GRF{i} ' left'],'Interpreter','none')
    xlabel('Gait cycle (%)')
end

subplot(3,3,4)
ylabel('GRF (% BW)')


subplot(3,3,[7:9])
hold on
axis equal
COP_r_x = R.ground_reaction.COP_r(R.ground_reaction.idx_stance_r,1);
COP_r_z = -R.ground_reaction.COP_r(R.ground_reaction.idx_stance_r,3);
plot(COP_r_x,COP_r_z,'.','color',colr,'DisplayName',legName)

COP_l_x = R.ground_reaction.COP_l(R.ground_reaction.idx_stance_l,1);
COP_l_z = -R.ground_reaction.COP_l(R.ground_reaction.idx_stance_l,3);
plot(COP_l_x,COP_l_z,'.','color',colr,'DisplayName',legName)

title('Centre of pressure in ground reference (top view)')
xlabel('x (m)')
ylabel('-z (m)')



