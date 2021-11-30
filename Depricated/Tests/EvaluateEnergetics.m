%% Evaluate energetics of our simulations
%----------------------------------------
% 

% file information
datapath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\Results\Simulation_Koelewijn2019';
% filename = 'Flat_s1Pog_speed_13_pp.mat';
% filename = 'Decl8_s1Pog_speed_13_pp.mat';
filename = 'Incl8_s1Pog_speed_13_pp.mat';

load(fullfile(datapath,filename),'R');

% Note this is a quit version: we should actuatlly do this on the
% non-interpolated data (at mesh and collocation points). in the function
% f_LoadSim_Gait92

% compute joint power
tau = R.Tid;
qd = R.Qdots;
Power = tau.*qd*pi/180;
TotalPower = sum(Power(:,7:end),2);

% compute the net work done during a gait cycle
Work = trapz(R.t,Power);
TotalWork = sum(Work(7:end));
disp(['Total work done: ' num2str(TotalWork) ' J']);
disp(['Total work done: ' num2str(TotalWork./R.StrideLength./R.body_mass) ' J/m/kg']);

% compute muscle work done
vM = -R.Muscle.vM;
Fce = R.Muscle.Fce;
MusclePower = Fce.*vM;
MuscleWorkInd = trapz(R.t,MusclePower);
MuscleWork = sum(MuscleWorkInd);
disp(['Total muscle work: ' num2str(MuscleWork) ' J']);
disp(['Total muscle work: ' num2str(MuscleWork./R.StrideLength./R.body_mass) ' J/m/kg']);

% work done by torque actuators
iTorqAct = [24:31];
Power_TorqueAct = tau(:,iTorqAct) .* qd(:,iTorqAct);
Work_TorqueActInd = trapz(R.t,Power_TorqueAct);
Work_TorqueAct = sum(Work_TorqueActInd);
disp(['Total work torque actuators: ' num2str(Work_TorqueAct) ' J']);
disp(['Total work torque actuators: ' num2str(Work_TorqueAct./R.StrideLength./R.body_mass) ' J/m/kg']);
