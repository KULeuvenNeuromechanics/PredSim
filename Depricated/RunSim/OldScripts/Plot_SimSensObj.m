%% Plot Simulations: 17/03/2020
%------------------------------
clear all; close all; clc;

%% Path information
Datapath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\Results';
DataFolders = {'Batch_SensObjective'};

S.LoadRes = 1;
S.Plot1   = 0;

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
if S.Plot1
    % visualise indexes with DispHeader
    DispHeader(Names);    
       
    % Plot one figure with all the lines
    ColV = [1 0 0; 0 1 0; 0 0 1];
    for i=1:nF
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

%% Analyse the specific results

DispHeader(Names);

% plot influence Obj A
IndV = 1:5;
CsV = copper(length(IndV));
ct = 1;
h = figure();
for i = IndV
    [path,name,ext] = fileparts(Names{i});
    PlotResults_3DSim(Names{i},CsV(ct,:),name,h); ct = ct+1;
end


% plot influence Obj A - Passive
IndV = 42:46;
CsV = copper(length(IndV));
ct = 1;
h = figure();
for i = IndV
    [path,name,ext] = fileparts(Names{i});
    PlotResults_3DSim(Names{i},CsV(ct,:),name,h); ct = ct+1;
end

% plot influence Obj E - active
IndV = 6:10;
CsV = copper(length(IndV));
ct = 1;
h = figure();
for i = IndV
    [path,name,ext] = fileparts(Names{i});
    PlotResults_3DSim(Names{i},CsV(ct,:),name,h); ct = ct+1;
end

% plot influence Obj E - NoExo
IndV = 27:31;
CsV = copper(length(IndV));
ct = 1;
h = figure();
for i = IndV
    [path,name,ext] = fileparts(Names{i});
    PlotResults_3DSim(Names{i},CsV(ct,:),name,h); ct = ct+1;
end

% plot influence Obj E - NoExo
IndV = 27:31;
CsV = copper(length(IndV));
ct = 1;
h = figure();
for i = IndV
    [path,name,ext] = fileparts(Names{i});
    PlotResults_3DSim(Names{i},CsV(ct,:),name,h); ct = ct+1;
end

% plot influence Obj Qdd assistance
IndV = 16:20;
CsV = copper(length(IndV));
ct = 1;
h = figure();
for i = IndV
    [path,name,ext] = fileparts(Names{i});
    PlotResults_3DSim(Names{i},CsV(ct,:),name,h); ct = ct+1;
end








