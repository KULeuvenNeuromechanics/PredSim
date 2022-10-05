
clear
clc
import casadi.*

% add repository to workspace
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
[pathRepoFolder,~,~] = fileparts(pathRepo);
results_folder = fullfile(pathRepoFolder,'PredSimResults');
addpath([pathRepo '/VariousFunctions'])

%%
result_path = fullfile([results_folder '\subject1_2D\subject1_2D_v10']);

load([result_path '.mat'],'R','model_info')

name_1 = 'Subject1_2D';
F  = external('F',fullfile(pathRepo,'Subjects',name_1,['F_' name_1 '.dll'])); 

data_ID = importdata(fullfile(pathRepo,'Subjects',name_1,'inverse_dynamics_v10.sto'));
data_IDa = importdata(fullfile(pathRepo,'Subjects',name_1,'inverse_dynamics_v10a.sto'));
data_IDb = importdata(fullfile(pathRepo,'Subjects',name_1,'inverse_dynamics_v10b.sto'));
data_IDc = importdata(fullfile(pathRepo,'Subjects',name_1,'inverse_dynamics_v10c.sto'));

colh = data_ID.colheaders(2:end)';
coord_names = model_info.ExtFunIO.coord_names.all;
idx1 = 1:10;
idx2 = [1:4,7,9,5,8,10,6];


T_ID = data_ID.data(1:100,2:end);
T_ID(:,idx1) = T_ID(:,idx2);

T_IDa = data_IDa.data(1:100,2:end);
T_IDa(:,idx1) = T_IDa(:,idx2);

T_IDb = data_IDa.data(1:100,2:end);
T_IDb(:,idx1) = T_IDb(:,idx2);

T_IDc = data_IDa.data(1:100,2:end);
T_IDc(:,idx1) = T_IDc(:,idx2);

qmot = read_motionFile_v40([result_path '.mot']);
Qs = qmot.data(1:100,2:11);

qmota = read_motionFile_v40([result_path 'a.mot']);
Qsa = qmota.data(1:100,2:11);

qmotb = read_motionFile_v40([result_path 'b.mot']);
Qsb = qmotb.data(1:100,2:11);

qmotc = read_motionFile_v40([result_path 'c.mot']);
Qsc = qmotc.data(1:100,2:11);

diff_q = R.kinematics.Qs - Qs;


diff1 = R.kinetics.T_ID - T_ID;
diff1b = R.kinetics.T_ID - T_IDb;

diff0 = R.kinetics.T_ID_0 - T_ID;

diff_T = R.kinetics.T_ID - R.kinetics.T_ID_0;

diff_ID = T_ID - T_IDc;

diff_IDb = T_IDb - T_IDa;

