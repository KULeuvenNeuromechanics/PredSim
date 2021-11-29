% get muscle parameters
%------------------------

clear all; close all; clc;
ModelPath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\OpenSimModel';
OsimModel = fullfile(ModelPath,'Rajagopal2015.osim');
Mnames = DispMusclesOsimModel(OsimModel);
muscleNames = Mnames(1:40);
MTparameters = ReadMuscleParameters(OsimModel,muscleNames);
save('MTparameters_Rajagopal2015.mat','MTparameters');
save('MuscleNames.mat','muscleNames');


%% save also an adapted version of the rajagopal model
% (see discussion in document (RajagopalModel.md/html)
iRF = find(strcmp(Mnames,'recfem_r'));
iVasInt = find(strcmp(Mnames,'vasint_r'));
iVasLat= find(strcmp(Mnames,'vaslat_r'));
iVasMed = find(strcmp(Mnames,'vasmed_r'));

MTparameters(2,iRF) = 0.1;
MTparameters(2,iVasInt) = 0.11;
MTparameters(2,iVasLat) = 0.11;
MTparameters(2,iVasMed) = 0.11;
save('MTparameters_Rajagopal2015_LongKneeM.mat','MTparameters');
save('MuscleNames.mat','muscleNames');
