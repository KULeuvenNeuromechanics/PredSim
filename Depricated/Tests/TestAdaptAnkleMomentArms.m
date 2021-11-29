%% Test adapting the moment arms for the ankle muscles
%------------------------------------------------------

% path information
MainPath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim';
S.PolyFolder = 's1_Poggensee';

% path with original casadifunctions
S.CasadiFunc_OriginFolders = 'Casadi_s1Pog_mtp'; 

% path with adapted moment arms
S.CasadiFunc_AdaptedFolders = 's1_Poggensee_dMAnkle1'; 

% import casadi libraries
import casadi.*

% load the casadifunctions
f_Orgin = Function.load(fullfile(MainPath,'CasADiFunctions',...
    S.CasadiFunc_OriginFolders,'f_lMT_vMT_dM'));

f_Adapt = Function.load(fullfile(MainPath,'CasADiFunctions',...
    S.CasadiFunc_AdaptedFolders,'f_lMT_vMT_dM'));

% sample input
q = zeros(1,10);
q(5) = 0.2;
qdot = zeros(1,10);

% evaluate polynomials
[LMT_or,~,dM_or] = f_Orgin(q,qdot);
[LMT_ad,~,dM_ad] = f_Adapt(q,qdot);

LMT = LMT_or-LMT_ad;
dM = dM_or-dM_ad;

% plot figure
figure();
subplot(1,2,1)
bar(full(LMT));
subplot(1,2,2)
bar(full(dM));