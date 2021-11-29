function [Tracking] = GetTrackingData(S,pathRepo,joints,N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

pathData = [pathRepo,'/OpenSimModel/',S.Tracking.subject];
trials =  S.Tracking.trials;
p = 1;
mtp_jointType = 'mtpPin';

nametrial(p).id    = [trials{p},'_',mtp_jointType]; % Experimental walking trial to track
nametrial(p).ID    = ['ID_',nametrial(p).id];
nametrial(p).GRF   = ['GRF_',nametrial(p).id];
nametrial(p).IK    = ['IK_',nametrial(p).id];

time_opt = load(fullfile(pathData,'IK',['Time_' trials{p} '.mat']));

% Extract joint kinematics
pathIK = [pathData,'/IK/',nametrial(p).IK,'.mat'];
Qs(p).p = getIK(pathIK,joints);
% Extract ground reaction forces and moments
pathGRF = [pathData,'/GRF/',nametrial(p).GRF,'.mat'];
GRF(p).p = getGRF(pathGRF);
% Extract joint kinetics
pathID = [pathData,'/ID/',nametrial(p).ID,'.mat'];
ID(p).p = getID(pathID,joints);
% Interpolation experimental data
t0 = time_opt.time(1);
tend = time_opt.time(2);
time_expi.ID(1) = find(round(ID(p).p.time,4) == t0);
time_expi.ID(2) = find(round(ID(p).p.time,4) == tend);
time_expi.GRF(1) = find(round(GRF(p).p.time,4) == t0);
time_expi.GRF(2) = find(round(GRF(p).p.time,4) == tend);
step = (ID(p).p.time(time_expi.ID(2))-ID(p).p.time(time_expi.ID(1)))/(N-1);
interval = t0:step:tend;
ID(p).p.allinterp = interp1(ID(p).p.all(:,1),ID(p).p.all,interval);
Qs(p).p.allinterpfilt = interp1(Qs(p).p.allfilt(:,1),Qs(p).p.allfilt,interval);
GRF(p).p.val.allinterp = interp1(round(GRF(p).p.val.all(:,1),4),...
    GRF(p).p.val.all,round(interval,4));
GRF(p).p.MorGF.allinterp = interp1(round(GRF(p).p.MorGF.all(:,1),4),...
    GRF(p).p.MorGF.all,round(interval,4));
GRF(p).p.pos.allinterp = interp1(round(GRF(p).p.pos.all(:,1),4),...
    GRF(p).p.pos.all,round(interval,4));
GRF(p).p.Mcop.allinterp = interp1(round(GRF(p).p.Mcop.all(:,1),4),...
    GRF(p).p.Mcop.all,round(interval,4));

% save in structure

Tracking.ID = ID;
Tracking.Qs = Qs;
Tracking.GRF = GRF;
Tracking.t = tend-t0;
Tracking.tWindow = time_opt.time;


end

