% Script to test the function getBounds
%
% Author: Lars D'Hondt
% Modified: Dhruv Gupta
% Date: 2 dec 2021
%
% TO DO
%--------------------------------------------------------------------------
clear
clc
close all

% add repository to workspace
% [pathHere,~,~] = fileparts(mfilename('fullpath'));
% [pathRepo,~,~] = fileparts(pathHere);
osim_path = 'test_1.osim';
% % load model info
% load('model_info.mat')

%% make tmp model_info
load('F_test_1_IO.mat');

model_info.ExtFunIO = IO;

import org.opensim.modeling.*;
model = Model(osim_path);
for i=1:model.getMuscles().getSize()
    muscle_names{i} = char(model.getMuscles().get(i-1).getName());
end

model_info.muscle_info.muscle_names = muscle_names;

%% create motion file
% colheaders
coordinate_names = fieldnames(model_info.ExtFunIO.coordi);
% coordinate_names_random = coordinate_names(randperm(length(coordinate_names)));
% q_in.labels = [{'time'},coordinate_names_random'];
q_in.labels = [{'time'},coordinate_names'];
% random data
q_in.data = [linspace(0.01,0.05,10)', (rand([10,length(coordinate_names)])-0.3)];
% angles are in radians
q_in.inDeg = 'no';
% filename
filename = 'test_IK.mot';
% generate motion file
write_motionFile_v40(q_in, filename)

%% create needed settings
% S.subject.folder_name = pwd;
S.subject.IKfile_bounds = 'test_IK.mot';
S.subject.vPelvis_x_trgt = 1.25;
S.subject.IG_pelvis_y = 0.85;
S.bounds.t_final.lower = 0.2;
S.bounds.t_final.upper = 2;
S.subject.mtp_type = 'active';
S.misc.gaitmotion_type = 'FullGaitCycle';

%%
Qs = getIK(S.subject.IKfile_bounds,model_info);

%% call getBounds
residualsi = 1:length(fields(model_info.ExtFunIO.coordi));
ground_pelvisi      = model_info.ExtFunIO.jointi.ground_pelvis; % ground-pelvis
trunki              = model_info.ExtFunIO.jointi.torso; % trunk
legsi               = [model_info.ExtFunIO.jointi.leg_r model_info.ExtFunIO.jointi.leg_l]; % arms
armsi               = [model_info.ExtFunIO.jointi.arm_r model_info.ExtFunIO.jointi.arm_l]; % arms
mtpi                = [model_info.ExtFunIO.jointi.mtp_r model_info.ExtFunIO.jointi.mtp_l]; % mtps
coord_noarmsi       = [ground_pelvisi legsi trunki]; % all but arms
coord_muscleActuated= [legsi trunki]; % all but arms
% Number of degrees of freedom for later use
nq.all      = length(residualsi); % all
nq.abs      = length(ground_pelvisi); % ground-pelvis
nq.trunk    = length(trunki); % trunk
nq.arms     = length(armsi); % arms
nq.mtp      = length(mtpi); % arms
nq.leg      = length(legsi);
nq.muscAct  = length(coord_muscleActuated);

model_info.ExtFunIO.nq = nq;
jointi = getJointi();
[bounds,scaling,Qs_spline] = getBounds_all(Qs,model_info,S);
[bounds_2,scaling_2] = getBounds_all_mtp(Qs,92,model_info.ExtFunIO.nq,jointi,S.subject.vPelvis_x_trgt);

time_IC         = [Qs.time(1),Qs.time(end)];
N=10;
d=3;
NMuscle=92;
guess = getGuess_DI_opti(Qs,N,time_IC,scaling,S,d,model_info);
% guess = getGuess_QR_opti(N,scaling,model_info,S,d);

visualizebounds



% %% NOT A GOOD TEST BECAUSE OF THE SPLINE FUNCTION
% for i=1:length(coordinate_names)
%     coord_name_i = coordinate_names{i};
%     q = q_in.data(:,find(strcmp(q_in.labels,coord_name_i)));
%     lb = bounds.Qs.lower(i);
%     ub = bounds.Qs.upper(i);
%     sf = scaling.Qs(i);
%     maxq = max(q);
%     minq = min(q);
%     if lb<(minq/sf) && (maxq/sf)<ub
%         check = 'ok';
%     else
%         check = 'nok';
%     end
% %     disp([coord_name_i '  ' check])
% 
%     disp([coord_name_i '    ' num2str(lb) ' < ' num2str(minq/sf) ' < ' num2str(maxq/sf) ' < ' num2str(ub) '  ' check])
% end

% for i=1:length(coordinate_names)
%     coord_name_i = coordinate_names{i};
%     q = q_in.data(:,strcmp(q_in.labels,coord_name_i));
%     lb = bounds_2.Qs.lower(i);
%     ub = bounds_2.Qs.upper(i);
%     sf = scaling_2.Qs(i);
%     maxq = max(q);
%     minq = min(q);
%     if lb<(minq/sf) && (maxq/sf)<ub
%         check = 'ok';
%     else
%         check = 'nok';
%     end
% %     disp([coord_name_i '  ' check])
% 
%     disp([coord_name_i '    ' num2str(lb) ' < ' num2str(minq/sf) ' < ' num2str(maxq/sf) ' < ' num2str(ub) '  ' check])
% end

%%
for i=1:length(coordinate_names)
    coord_name_i = coordinate_names{i};
    lb = bounds.Qs.lower(i);
    ub = bounds.Qs.upper(i);
    lb_2 = bounds_2.Qs.lower(i);
    ub_2 = bounds_2.Qs.upper(i);
    sf = scaling.Qs(i);
    sf_2 = scaling_2.Qs(i);
    
    lbDiff = lb-lb_2;
    ubDiff = ub-ub_2;
    sfDiff = sf-sf_2;
    if lb==lb_2 && ub_2==ub && sf_2==sf
        check = 'ok';
    else
        check = 'nok';
    end
    disp([coord_name_i '    ' num2str(lbDiff) ' ' num2str(ubDiff) ' ' num2str(sfDiff) ' ' check])
    
%     disp([coord_name_i '    ' num2str(lb) ' < ' num2str(minq/sf) ' < ' num2str(maxq/sf) ' < ' num2str(ub) '  ' check])
    
end
