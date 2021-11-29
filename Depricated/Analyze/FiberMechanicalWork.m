
%% Detailed analysis of metabolic energy
%----------------------------------------

% datapath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\Results\ExpTorques_Lower_lba001';
% ResultsFile = 'NoExo.mat';


[R] = f_LoadSim_Gait92('ExpTorques_Lower_lba001','NoExo');

%% Test if computed muscle fiber lengths make sense

nfr = length(R.t);
lMtilde = R.lMtilde;
lMo = R.Muscle.MTparameters(2,:);
lMo_m = ones(nfr,1) * lMo;
lM = lMtilde .* lMo_m;

lM_dot = diff(lM)./diff(R.t);

iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));
figure(); 
subplot(1,2,1)
plot(R.t,lM(:,iSol));
xlabel('Time [s]');
ylabel('Fiber length');

subplot(1,2,2)
plot(R.t,R.Muscle.vM(:,iSol)); hold on;
plot(R.t(1:end-1),lM_dot(:,iSol),'--r'); hold on;

xlabel('Time [s]');
ylabel('Fiber velocity');

%% Compute positive mechanical work of muscle fibers

MusclePower = -R.Muscle.vM .* R.Muscle.Fce;

% plot fiber power of a selected muscles
P_msel = MusclePower(:,strcmp(R.colheaders.muscles,'lat_gas_r'));
figure();
plot(P_msel);

% find muscles with high peak power
Maxpower = max(MusclePower);
iSel = Maxpower>10;
DispHeader(R.colheaders.muscles(iSel))
figure();
plot(MusclePower(:,iSel))

% get muscle work
FiberWork = trapz(R.t,MusclePower);
sum(FiberWork);

% get positive muscle work
MusclePowerPos = MusclePower;
MusclePowerPos(MusclePowerPos<0) = 0;
PosFiberWork = trapz(R.t,MusclePowerPos);
PosFiberWork_total = sum(PosFiberWork);
% average power
PosPower = PosFiberWork_total./(R.t(end)-R.t(1));
disp(['Average positive power during walking: ' num2str(PosPower) ' W']);

% get negative muscle work
MusclePowerNeg = MusclePower;
MusclePowerNeg(MusclePowerNeg>0) = 0;
NegFiberWork = trapz(R.t,MusclePowerNeg);
NegFiberWork_total = sum(NegFiberWork);
% average power
NegPower = NegFiberWork_total./(R.t(end)-R.t(1));
disp(['Average negative power during walking: ' num2str(NegPower) ' W']);

% simple estimate of metabolic work:
SimpleMetabE = NegPower / -1.2 + PosPower / 0.25;
disp(['Simple estimate metabolic power: ' num2str(SimpleMetabE), ' W']);

%% average metabolic power of different models

Names = fieldnames(R.Energy);
for i=1:length(Names)
    Psel = nanmean(R.Energy.(Names{i}));
    disp(['P metab avg ' Names{i} '  ' num2str(Psel) ' W']);    
end



