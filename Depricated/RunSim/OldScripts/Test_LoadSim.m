% TestLoadSim function

ResultsFolder = 'TestS1_PoggenSee';
OutName = 'Assistance_100';
R = f_LoadSim_PoggenSee2020_DefaultS(ResultsFolder,OutName);


ResultsFolder = 'TestS1_PoggenSee';
OutName = 'NoExo_Stiffmtp';
f_LoadSim_PoggenSee2020_DefaultS(ResultsFolder,OutName);

% ResultsFolder = 'TestS1_PoggenSee';
% OutName = 'NoExo_NormalMTP_pp';
% f_LoadSim_PoggenSee2020_DefaultS(ResultsFolder,OutName);

ResultsFolder = 'TestS1_PoggenSee';
OutName = 'NoExo';
f_LoadSim_PoggenSee2020_DefaultS(ResultsFolder,OutName);

% plot results

%% Plot results

Rpath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\Results';
ResultsFolder = 'TestS1_PoggenSee';

Rname = fullfile(Rpath,ResultsFolder,'NoExo_Stiffmtp_pp.mat');
PlotResults_3DSim(Rname,[0 0 0]);
Rname = fullfile(Rpath,ResultsFolder,'NoAssistance_pp.mat');
PlotResults_3DSim(Rname,[0 0 1],gcf);
Rname = fullfile(Rpath,ResultsFolder,'Assistance_100_pp.mat');
PlotResults_3DSim(Rname,[1 0 0],gcf);

