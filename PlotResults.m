%% Plot solution


% plot results
ResPath = 'C:\Users\mat950\Documents\Software\Sim\PredSim\Results\Falisse_et_al_2022';
PathExpData = 'C:\Users\mat950\Documents\Software\Sim\PredInt\src\ExpData\ExperimentalData.mat';
h = figure;
ResultsFile = fullfile(ResPath,'Falisse_et_al_2022_v3.mat');
% get the color
Cs = [0, 0, 1];
% legend name
LegName = 'SimExo';
% x-axis info
xVal = 1;
xLab = 'Sim';
ExpSubj = {'subject1'};
% plot on the current figure
PlotResults_3DSim(ResultsFile,Cs,LegName,PathExpData,h,xVal,xLab,ExpSubj);

figure();
load(ResultsFile,'R');
plot(R.time.mesh_GC(:,1:end-1), R.exo.Texo);
