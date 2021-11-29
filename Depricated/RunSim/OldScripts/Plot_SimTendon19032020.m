%% Plot Simulations: 17/03/2020
%------------------------------
clear all; close all; clc;

%% Path information
Datapath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\Results';
DataFolders = {'Batch_TendonStiff'};

S.LoadRes = 1;
S.Plot    = 0;

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
    ColV1 = copper(5);
    ColV2 = jet(5);
    h1 = figure();
    for i=1:5
        FileSel = Names{i};
        [path,name,ext] = fileparts(FileSel);
        PlotResults_3DSim(FileSel,ColV1(i,:),name,h1);
        FileSel = Names{i+6};
        [path,name,ext] = fileparts(FileSel);
        PlotResults_3DSim(FileSel,ColV2(i,:),name,h1);
    end
    
    % plot normal walking with exo assistance
    TendonK = [15 20 25 30 35];
    h2 = figure();
    for i=1:5
        FileSel = Names{i};
        [path,name,ext] = fileparts(FileSel);
        PlotResults_3DSim(FileSel,ColV1(i,:),name,h2,TendonK(i),'Stiffness');
    end
    
    
    % plot normal walking withou exo assistance
    TendonK = [15 20 25 30 35];
    h3 = figure();
    for i=1:5
        FileSel = Names{i+6};
        [path,name,ext] = fileparts(FileSel);
        PlotResults_3DSim(FileSel,ColV1(i,:),name,h3,TendonK(i),'Stiffness');
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

