import casadi.*

%% Path information
pathmain = pwd;
% We use different external functions, since we also want to access some
% parameters of the model in a post-processing phase.
[pathRepo,~,~] = fileparts(pathmain);
pathExternalFunctions = [pathRepo,'/ExternalFunctions'];
pathExtSource = fullfile(pathExternalFunctions,'GitSource');
% Loading external functions.

%% Copy files from opensim force

Files = {'SimExo_3D.ilk','SimExo_3D.dll','SimExo_3D.pdb','SimExo_3D.exp','SimExo_3D.lib'};
OsimBuild = 'C:\OpenSimGit\opensim-core-build\RelWithDebInfo';
for i=1:length(Files)
    copyfile(fullfile(OsimBuild,Files{i}),fullfile(pathExtSource,Files{i}));
end

cd(pathExternalFunctions);
F1  = external('F','SimExo_3D_talus.dll');
F3  = external('F','SimExo_3D_ExportAll.dll');
cd(OsimBuild);
F2 = external('F','SimExo_3D.dll',struct('enable_fd',true,...
    'enable_forward',false,'enable_reverse',false,...
    'enable_jacobian',false,'fd_method','backward'));
cd(pathmain);

%% load a decent state from a previous solution (with an old .dll file)


load('C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\Results\Sens_StepWidth\Passive_dCalcn_14cm_pp.mat');

nfr = length(R.Qs(:,1));
qdqdd = zeros(nfr,62);
qdqdd(:,1:2:62) = R.Qs;
qdqdd(:,2:2:62) = R.Qdots;
qdd = R.Qddots;
qdqdd(:,[1:6 13:end]) = qdqdd(:,[1:6 13:end])*pi./180;

qdqdd(:,11) = qdqdd(:,11);

COPR = zeros(nfr,3);
COPL = zeros(nfr,3);
FR = zeros(nfr,3);
FL = zeros(nfr,3);
MR = zeros(nfr,3);
ML = zeros(nfr,3);
Out_F2Vect = zeros(nfr,85);

for ind = 1:nfr
    Out_F1 = full(F1([qdqdd(ind,:)'; qdd(ind,:)'; 0;0]));
    Out_F2 = full(F2(qdqdd(ind,:)',[qdd(ind,:)'; 0;0]));    
    Out_F2Vect(ind,:) = Out_F2;
    % compute the COP position    
    FR(ind,:) = -Out_F2(68:70);
    FL(ind,:) = -Out_F2(71:73);
    MR(ind,:) = -Out_F2(74:76);
    ML(ind,:) = -Out_F2(77:79);
    if abs(FR(ind,2)) > 10
        COPR(ind,:) = [MR(ind,3)./FR(ind,2), 0, -MR(ind,1)./FR(ind,2)];
    end
    if abs(FL(ind,2)) > 10
        COPL(ind,:) = [ML(ind,3)./FL(ind,2), 0, -ML(ind,1)./FL(ind,2)];
    end
end

figure();
subplot(1,2,1);
plot(COPL(:,1),'b'); hold on; plot(COPR(:,1),'r');
xlabel('frames'); ylabel('COPx [m]');
subplot(1,2,2);
plot(COPL(:,3),'b'); hold on; plot(COPR(:,3),'r');
xlabel('frames'); ylabel('COPz [m]');
legend('left foot','right foot');


%%
COPR2 = zeros(nfr,3);
COPL2 = zeros(nfr,3);
FR2 = zeros(nfr,3);
FL2 = zeros(nfr,3);
MR2 = zeros(nfr,3);
ML2 = zeros(nfr,3);
Out_F3Vect = zeros(nfr,73);
for ind = 1:nfr
    Out_F3 = full(F3([qdqdd(ind,:)'; qdd(ind,:)'; 0;0]));    
    Out_F3Vect(ind,:) = Out_F3;
    % compute the COP position    
    FR2(ind,:) = Out_F3(32:34);
    FL2(ind,:) = Out_F3(35:37);
    MR2(ind,:) = Out_F3(68:70);
    ML2(ind,:) = Out_F3(71:73);
    if abs(FR(ind,2)) > 10
        COPR2(ind,:) = [MR2(ind,3)./FR2(ind,2), 0, -MR2(ind,1)./FR2(ind,2)];
    end
    if abs(FL(ind,2)) > 10
        COPL2(ind,:) = [ML2(ind,3)./FL2(ind,2), 0, -ML2(ind,1)./FL2(ind,2)];
    end
end

subplot(1,2,1);
plot(COPL2(:,1),'--k');
plot(COPR2(:,1),'--k');
subplot(1,2,2);
plot(COPL2(:,3),'--k');
plot(COPR2(:,3),'--k');
legend('left foot','right foot');

% figure(); plot(FR2); hold on;
% plot(Out_F3Vect(:,32:34),'--k');

%% Export to mot file
% 
% data = [R.t FR COPR FL COPL zeros(nfr,6)];
% colnames = get_GRFlabels();
% filename =  'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\Results\Sens_StepWidth\Passive_dCalcn_14cm_pp_GRF.mot';
% generateMotFile(data, ['time ' colnames], filename);



% figure(); plot(FR); hold on;
% plot(Out_F2Vect(:,32:34),'--k');

% % figure(); plot(FR); hold on;
% figure();
% plot(Out_F2Vect(:,74:76),'b'); hold on;
% plot(Out_F2Vect(:,80:82),'--k');

%% Test computation wrench


figure();
subplot(2,2,1)
plot(FR); hold on; plot(FR2,'--k');
subplot(2,2,2)
plot(FL); hold on; plot(FL2,'--k');
subplot(2,2,3)
plot(MR); hold on; plot(MR2,'--k');
subplot(2,2,4)
plot(ML); hold on; plot(ML2,'--k');
