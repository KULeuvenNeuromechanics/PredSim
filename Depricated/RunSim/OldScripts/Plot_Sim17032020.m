%% Plot Simulations: 17/03/2020
%------------------------------
clear all; close all; clc;

%% Path information
Datapath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\Results';
DataFolders = {'BatchSim_2020_03_16_e','BatchSim_2020_03_16_UpdIG','BatchSim_2020_03_17_UpdIG'};

S.LoadRes = 1;
S.Plot    = 1;

%% Post process the simulations

ct = 1;
nF = length(DataFolders);
for f = 1:nF
    % get all the simulation results in this folder
    dpath = fullfile(Datapath,DataFolders{f});
    MatFiles = dir(fullfile(dpath,'*.mat'));
    nFil = length(MatFiles);
    for i = 1:nFil
        filename = MatFiles(i).name;
        FileEnd = filename(end-6:end);
        OutName = fullfile(dpath,[filename(1:end-4) '_pp.mat']);
        if ~strcmp(FileEnd,'_pp.mat')
            Names{ct} = OutName;
            FolderIndex(ct) = f;
            ct= ct+1;
            if S.LoadRes && ~exist(OutName,'file')
                f_LoadSim_PoggenSee2020_DefaultS(DataFolders{f},filename);
                disp(ct);
            end
        end
    end
end


%% Plot the simulation results for analysis
if S.Plot
    % visualise indexes with DispHeader
    DispHeader(Names);
    
    % plot for first folder
    dpath = fullfile(Datapath,DataFolders{1});
    PlotResults_3DSim(Names{18},[0 0 1],'Passive'); h = gcf;
    PlotResults_3DSim(Names{17},[0 0 0],'NoExo',h);
    Cols = copper(16);
    for i=1:16
        FileSel = Names{i};
        [path,name,ext] = fileparts(FileSel);
        PlotResults_3DSim(FileSel,Cols(i,:),name,h);
    end
    % saveas(h,fullfile(dpath,'FigureResults.fig'));
    
    % plot for second folder
    dpath = fullfile(Datapath,DataFolders{2});
    PlotResults_3DSim(Names{35},[0 0 1],'Passive'); h = gcf;
    PlotResults_3DSim(Names{34},[0 0 0],'NoExo',h);
    Inds = 19:33; Cols = copper(length(Inds)); ct = 1;
    for i=Inds
        FileSel = Names{i};
        [path,name,ext] = fileparts(FileSel);
        PlotResults_3DSim(FileSel,Cols(ct,:),name,h); ct = ct+1;
    end
    % saveas(h,fullfile(dpath,'FigureResults.fig'));
    
    
    % plot for third folder
    dpath = fullfile(Datapath,DataFolders{3});
    PlotResults_3DSim(Names{59},[0 0 1],'Passive'); h = gcf;
    PlotResults_3DSim(Names{58},[0 0 0],'NoExo',h);
    Inds = 36:57; Cols = copper(length(Inds)); ct = 1;
    for i=Inds
        FileSel = Names{i};
        [path,name,ext] = fileparts(FileSel);
        PlotResults_3DSim(FileSel,Cols(ct,:),name,h); ct = ct+1;
    end
    % saveas(h,fullfile(dpath,'FigureResults.fig'));
    
    % Plot one figure with all the lines
    ColV = [1 0 0; 0 1 0; 0 0 1];
    for i=1:3
        Inds = find(FolderIndex ==i);
        for j = Inds
            FileSel = Names{j};
            [path,name,ext] = fileparts(FileSel);
            if i == 1 && j==1
                PlotResults_3DSim(FileSel,ColV(i,:),name);
            else
                PlotResults_3DSim(FileSel,ColV(i,:),name,gcf);
            end
        end
    end
end

%% Plot the relevant figure for the report

% Inds = [46 47 34 19];
% HeadNam = {'Without','Passive','Active','Active'};
% ColV = [0 0 0; 1 0 0; 0 0 1; 0 0 1];
%
% for i=length(Inds)
%     FileSel = Names{Inds(i)};
%     if i== 1
%         PlotResults_3DSim(FileSel,ColV(i,:),HeadNam{i});
%     else
%         PlotResults_3DSim(FileSel,ColV(i,:),HeadNam{i},gcf);
%     end
% end

