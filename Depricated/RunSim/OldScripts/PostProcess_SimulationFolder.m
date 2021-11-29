%% PostProcess Simluations
%---------------------------


% Default script for post-processing of all simulations
%------------------------------------------------------
clear all; clc;


%% Path information
Datapath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\Results';
% DataFolders = {'BatchSim_2020_03_16_e','BatchSim_2020_03_16_UpdIG','BatchSim_2020_03_17_UpdIG',...
%     'Batch_TendonStiff','Batch_SensObjective'};
% DataFolders = {'BatchSim_2020_03_17_UpdIG'};
DataFolders = {'Batch_TendonScale_k15','Batch_TendonScale_k20','Batch_TendonScale_k25',...
    'Batch_TendonScale_k30','Batch_TendonScale_k35','Batch_TendonScale_k40','BatchSim_2020_03_17_UpdIG','Batch_TendonStiff','Batch_SensObjective'};

S.OverWrite = 1;

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
            if (~exist(OutName,'file') || S.OverWrite == 1)
                f_LoadSim_PoggenSee2020_DefaultS(DataFolders{f},filename);
%                 disp(ct);
            end
        end
    end
end