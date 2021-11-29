%% Test if we read the MT parameters correctly


Mainpath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim';

% test MT properties
New = load(fullfile(Mainpath,'CasADiFunctions\s1_Poggensee_FisoRaja\MTparameters.mat'));
Old = load(fullfile(Mainpath, 'MuscleModel\s1_Poggensee\MTparameters_s1_Poggensee_mtp.mat'));


DiffMT = New.MTparameters(3,1:46)-Old.MTparameters(3,1:46);
figure(); bar(DiffMT);

%% Test casadi functions

% Fpath1 = 's1_Poggensee_FisoRaja';
Fpath1 = 's1_Poggensee_testLars';
Fpath2 = 'Casadi_s1Pog_mtp';


% test muscle tendon length
import casadi.*
CasPath1 = fullfile(Mainpath,'CasADiFunctions',Fpath1);
f_lMT_vMT_dM_New = Function.load(fullfile(CasPath1,'f_lMT_vMT_dM'));
CasPath2 = fullfile(Mainpath,'CasADiFunctions',Fpath2);
f_lMT_vMT_dM_Old = Function.load(fullfile(CasPath2,'f_lMT_vMT_dM'));

q = -rand(1,10);
qd = -rand(1,10);

[LMT1,vMT1,dM1] = f_lMT_vMT_dM_New(q,qd);
[LMT2,vMT2,dM2] = f_lMT_vMT_dM_Old(q,qd);

dLMT = full(LMT1)-full(LMT2);
dVMT = full(vMT1)-full(vMT2);
dM = full(dM1)-full(dM2);

max(abs(diff(dLMT)))
max(abs(diff(dVMT)))
max(abs(diff(dM)))

% test other functions ?


%% check mass file
% 
% Mdef = load(fullfile(Mainpath,'CasADiFunctions','MassM.mat'));
% Mcas = load(fullfile(Mainpath,'CasADiFunctions',Fpath1,'MassM.mat'));
% 
% disp(num2str(max(Mdef.MassM-Mcas.MassM)));