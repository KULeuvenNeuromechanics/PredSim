function [fig_hand] = plot_figure_E_muscle_bar(R,varargin)
% Original author: Lars D'Hondt
% Original date: 24/May/2022
%
% Last edit by: Tom Buurke
% Last edit date: 27/Oct/2023
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
    elseif isa(varargin{i},'struct')
        model_info = varargin{i};
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

%% Make plot E Muscle Bar
iM = find(contains(R.colheaders.muscles,'_r'));

mVect = R.colheaders.muscles(iM);

mus_cat = categorical(mVect);
Emus = zeros(length(mus_cat),1);

mus_cat = reordercats(mus_cat,[length(mVect):-1:1]);


imu_ct = 1;
for imu=1:length(mVect)
    Emus(imu) = trapz(R.time.mesh_GC(1:end-1),R.metabolics.Bhargava2004.Edot_gait(:,iM(imu_ct)))/model_info.mass/R.spatiotemp.dist_trav;
    imu_ct = imu_ct+1;
end

hold on
br1=barh(mus_cat,Emus);
grid on
for ibr=1:length(br1)
    br1(ibr).FaceColor = 'flat';
    br1(ibr).CData = colr;
    br1(ibr).FaceAlpha=0.5;
    br1(ibr).EdgeAlpha=0.5;
end
title('Metabolic energy of muscles')
xlabel('E_{metab} (Jkg^{-1}m^{-1})')
    lh24=legend('Nominal','Asymmetric','location','northeast','NumColumns',2);
    lhPos = lh24.Position;
    lhPos(1) = lhPos(1)+0.1;
    set(lh24,'position',lhPos);
end