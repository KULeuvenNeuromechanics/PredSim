import casadi.*

pathmain = pwd;
% We use different external functions, since we also want to access some
% parameters of the model in a post-processing phase.
[pathRepo,~,~] = fileparts(pathmain);
pathExternalFunctions = [pathRepo,'/ExternalFunctions'];
% Loading external functions.
cd(pathExternalFunctions);
setup.derivatives =  'AD'; % Algorithmic differentiation
F  = external('F','SimExo_3D_talus.dll');
Fs  = external('F','SimExo_3D_talus_small.dll');
F2  = external('F','SimExo_3D_talus_out.dll');
cd(pathmain);

%% load a decent state from a previous solution (with an old .dll file)


load('C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\Results\BatchSim_2020_03_17_UpdIG\Passive_pp.mat');

qdqdd = zeros(length(R.Qs(:,1)),62);
qdqdd(:,1:2:62) = R.Qs;
qdqdd(:,2:2:62) = R.Qdots;
qdd = R.Qddots;
ind = 15;
qdqdd(:,[1:6 13:end]) = qdqdd(:,[1:6 13:end])*pi./180;


Texo = [0 0]'; % left and right leg
Tpas = full(F([qdqdd(ind,:)'; qdd(ind,:)'; Texo]));
Texo = [10 10]'; % left and right leg
Tact10 = full(F([qdqdd(ind,:)'; qdd(ind,:)'; Texo]));

disp(Tpas');
disp(Tact10'-Tpas');

%% Evaluate if F and F2 give the same results

T1 = full(F([qdqdd(ind,:)'; qdd(ind,:)'; 10;10;]));
T2 = full(F2([qdqdd(ind,:)'; qdd(ind,:)'; 10;10;]));
figure(); plot(T1(1:31)-T2(1:31),'*k');
disp(sum(abs(T1(1:31)-T2(1:31))));

