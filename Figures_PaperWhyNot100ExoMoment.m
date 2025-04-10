%% Make figure 4 paper
RPath = 'C:\Users\mat950\OneDrive - Vrije Universiteit Amsterdam\Onderzoek\SimResults\Res_pcAnkleExo';

addpath(genpath(pwd));
SetFigureDefaults();
IDVect = 0:10:90;
h = figure();
% Cols = copper(10); % To Do: adapt this
[gradient_colors_red, gradient_colors_exo, gradient_colors_green, ...
    gradient_colors_eDot, gradient_colors_exp] = getcolorgrad_israel(10);
Cols = gradient_colors_eDot;
mk = 6;
Coords = {'ankle_angle_r','knee_angle_r','hip_flexion_r'};
title_coord = {'Ankle','Knee','Hip'};
ModelSim = 'Dhondt2023';
% ModelSim = 'Falisse_et_al_2022';
for i=1:length(IDVect)
    ResPath = fullfile(RPath,['Results_exo_ID_' num2str(IDVect(i))],ModelSim);
    MatFiles = dir(fullfile(ResPath,'*.mat'));
    ResFile = fullfile(ResPath,MatFiles(1).name);
    load(ResFile,'R','model_info','setup');
    Cs = Cols(i,:);
    for j = 1:length(Coords)
        subplot(3,3,j)
        plot(R.kinematics.Qs(:,strcmp(R.colheaders.coordinates,Coords{j})),'Color',Cs); hold on;
        if j ==1
            ylabel('joint angle [deg]')
        end
        title(title_coord{j});
    end
    for j = 1:length(Coords)
        subplot(3,3,j+3)
        plot(R.kinetics.T_ID(:,strcmp(R.colheaders.coordinates,Coords{j})),'Color',Cs); hold on;
        if j ==1
            ylabel('joint moment [Nm]')
        end
        title(title_coord{j});
    end

    % exoskeleton power
    subplot(3,3,7)
    Tankle = R.kinetics.T_ID(:,strcmp(R.colheaders.coordinates,'ankle_angle_r'));
    qdankle = R.kinematics.Qdots(:,strcmp(R.colheaders.coordinates,'ankle_angle_r'))*pi/180;
    if isfield(R,'exo')
        Texo = R.exo.tExo_gc(:,2);
    elseif isfield(R.S,'AnklePctID') && R.S.AnklePctID
        Texo = Tankle.* R.S.AnklePctID_value;
    else
        Texo = 0;
    end
    Pexo = Texo.*qdankle;
    plot(Pexo,'Color',Cs); hold on;
    ylabel('exoskeleton power [W]');
    title('exo power');

    % plot activation and metabolic cost
    subplot(3,3,8)
    act_sq = sumsqr(R.muscles.a)./(100);
    act = mean(sum(R.muscles.a,2)./92);
    yyaxis left
    Cs_mus = gradient_colors_red(i,:);
    plot(IDVect(i),act,'o','Color',Cs_mus,'MarkerFaceColor',Cs_mus,'MarkerSize',mk); hold on;
    ylabel('muscle activation');
    yyaxis right
    Cs_metab = gradient_colors_eDot(i,:);
    MetabPower = mean(sum(R.metabolics.Bhargava2004.Edot_gait,2)./92);
    plot(IDVect(i),MetabPower,'o','Color',Cs_metab,'MarkerFaceColor',Cs_metab,'MarkerSize',mk); hold on;
    ylabel('metabolic power');
    title('effort');

    % plot objective value
    subplot(3,3,9)
    plot(IDVect(i),sum(R.objective.absoluteValues),'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
    ylabel('objective value');
    title('objective');
end

for i =1:9
    subplot(3,3,i)
    set(gca,'box','off');
end

%% overview figure -- Falisse model

ResPath = 'C:\Users\mat950\Documents\Software\Sim\PredSim\Results';
PathExpData = 'C:\Users\mat950\Documents\Software\Sim\PredInt\src\ExpData\ExperimentalData.mat';

RPath = 'C:\Users\mat950\Documents\Software\Sim\PredSim\Results';
h = figure();
IDVect = 20:10:80;
Cols = copper(10);
for i=1:length(IDVect)
    ResPath = fullfile(RPath,['Results_exo_ID_' num2str(IDVect(i))],'Falisse_et_al_2022');
    MatFiles = dir(fullfile(ResPath,'*.mat'));
    ResFile = fullfile(ResPath,MatFiles(1).name);
    % set the color
    Cs = Cols(i,:);
    % legend name
    LegName = ['Exo_ID_' num2str(IDVect(i))];
    % x-axis info
    xVal = IDVect(i);
    xLab = 'Exo Assistance';
    ExpSubj = {'subject1'};
    % plot on the current figure
    PlotResults_3DSim(ResFile,Cs,LegName,PathExpData,h,xVal,xLab,ExpSubj);
end
% add without exo
ResFile = 'C:\Users\mat950\Documents\Software\Sim\PredSim\Results\Fal22Ref_125\Falisse_et_al_2022_job6613.mat';
% set the color
Cs =  [0.5 0.5 1];
% legend name
LegName = 'Normal';
% x-axis info
xVal = 0;
xLab = 'Exo Assistance';
ExpSubj = {'subject1'};
% plot on the current figure
PlotResults_3DSim(ResFile,Cs,LegName,PathExpData,h,xVal,xLab,ExpSubj);

%% overview figure -- Dhondt model

ResPath = 'C:\Users\mat950\Documents\Software\Sim\PredSim\Results';
PathExpData = 'C:\Users\mat950\Documents\Software\Sim\PredInt\src\ExpData\ExperimentalData.mat';
h = figure();
IDVect = 20:10:80;
Cols = copper(10);
for i=1:length(IDVect)
    ResPath = fullfile(RPath,['Results_exo_ID_' num2str(IDVect(i))],'Dhondt2023');
    MatFiles = dir(fullfile(ResPath,'*.mat'));
    ResFile = fullfile(ResPath,MatFiles(1).name);
    % set the color
    Cs = Cols(i,:);
    % legend name
    LegName = ['Exo_ID_' num2str(IDVect(i))];
    % x-axis info
    xVal = IDVect(i);
    xLab = 'Exo Assistance';
    ExpSubj = {'subject1'};
    % plot on the current figure
    PlotResults_3DSim(ResFile,Cs,LegName,PathExpData,h,xVal,xLab,ExpSubj);
end
% add without exo
ResFile = fullfile(RPath,'Results_exo_ID_0\Dhondt2023\Dhondt2023_job162.mat');
% set the color
Cs =  [0.5 0.5 1];
% legend name
LegName = 'Normal';
% x-axis info
xVal = 0;
xLab = 'Exo Assistance';
ExpSubj = {'subject1'};
% plot on the current figure
PlotResults_3DSim(ResFile,Cs,LegName,PathExpData,h,xVal,xLab,ExpSubj);

% % add optimal assistance to this graph (to check if this is really the optimal solution)
% % add without exo
% ResFile = 'C:\Users\mat950\Documents\Software\Sim\PredSim\Results\Results_exo\Falisse_et_al_2022\Falisse_et_al_2022_job130.mat';
% % set the color
% Cs =  [1 0.5 0.5];
% % legend name
% LegName = 'OptimalExo';
% % x-axis info
% xVal = 0;
% xLab = 'Exo Assistance';
% ExpSubj = {'subject1'};
% % plot on the current figure
% PlotResults_3DSim(ResFile,Cs,LegName,PathExpData,h,xVal,xLab,ExpSubj);