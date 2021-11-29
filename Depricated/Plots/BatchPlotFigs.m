

dpath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\Results';
OutPath = 'C:\Users\u0088756\Box Sync\SimExo\DefaultFigures';

FolderNames = {'MainResults','BatchSim_2020_03_17_UpdIG','Batch_SensObjective','Batch_TendonStiff',...
    'Batch_TendonScale_k15','Batch_TendonScale_k20','Batch_TendonScale_k25','Batch_TendonScale_k30',...
    'Batch_TendonScale_k35','Batch_TendonScale_k40'};

OutFNames = {'MainResults','Sensitivity_ScalingExoTorque','Sensitivity_Objective','Sensitivity_TendonStiff',...
    'Sensitivity_TendonStiff_ExoScale/k15','Sensitivity_TendonStiff_ExoScale/k20','Sensitivity_TendonStiff_ExoScale/k25','Sensitivity_TendonStiff_ExoScale/k30',...
    'Sensitivity_TendonStiff_ExoScale/k35','Sensitivity_TendonStiff_ExoScale/k40'};

nf= length(FolderNames);

for i= 1:nf
    Fsel = fullfile(dpath,FolderNames{i});
    Fout = fullfile(OutPath,OutFNames{i});
    Plot3D_pwd(Fsel);
    OutFile = fullfile(Fout,'FigureResults.fig');
    if exist(OutFile,'file')
        delete(OutFile);
    end
    copyfile(fullfile(Fsel,'FigureResults.fig'),fullfile(Fout,'FigureResults.fig'));
end