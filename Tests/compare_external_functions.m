
% Lars D'Hondt
% 17/May/2022


clear
clc
import casadi.*

% add repository to workspace
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/PreProcessing'])

% load 1st 
name_1 = 'Subject1_2D';
F1  = external('F',fullfile(pathRepo,'Subjects',name_1,['F_' name_1 '.dll'])); 
JF1  = external('jac_F',fullfile(pathRepo,'Subjects',name_1,['F_' name_1 '.dll'])); 

load(fullfile(pathRepo,'Subjects',name_1,['F_' name_1 '_IO.mat'])); 
IO1 = IO;



% load 2nd
name_2 = 'Subject1_2D_3D';
F2  = external('F',fullfile(pathRepo,'Subjects',name_2,['F_' name_2 '.dll'])); 
JF2  = external('jac_F',fullfile(pathRepo,'Subjects',name_2,['F_' name_2 '.dll'])); 

load(fullfile(pathRepo,'Subjects',name_2,['F_' name_2 '_IO.mat'])); 
IO2 = IO;

% F2 = external('F','C:\Users\u0150099\OneDrive - KU Leuven\PhD\algodiff\predictiveSimulations_2D\ExternalFunctions\PredSim_2D_pp.dll');
% IO2.coordi.pelvis_tilt  = 1; 
% IO2.coordi.pelvis_tx    = 2;
% IO2.coordi.pelvis_ty    = 3;
% IO2.coordi.hip_flexion_l        = 4;
% IO2.coordi.hip_flexion_r        = 5;
% IO2.coordi.knee_angle_l       = 6;
% IO2.coordi.knee_angle_r       = 7;
% IO2.coordi.ankle_angle_l      = 8;
% IO2.coordi.ankle_angle_r      = 9;
% IO2.coordi.lumbar_extension    = 10;
% IO2.GRFs.right_foot = [11,12];
% IO2.GRFs.left_foot = [13,14];


%%


coord_names = fieldnames(IO1.coordi);
n_coord = length(coord_names);
n_coord2 = length(fieldnames(IO2.coordi));


idx2 = zeros(n_coord,1);
idx1 = (1:n_coord)';
for i=1:n_coord
    idx2(i) = IO2.coordi.(coord_names{i});
end

% coord_names_all = [coord_names,fieldnames(IO2.coordi)];

idx1_r = [];
idx1_l = [];

for i=1:n_coord
    cni = coord_names{i};
    if strcmp(cni(end-1:end),'_r')
        idx1_r(end+1) = i;
        idx1_l(end+1) = IO1.coordi.([cni(1:end-1) 'l']);
    end

end

IO1_idx = [idx1',IO1.GRFs.right_foot(1:2),IO1.GRFs.left_foot(1:2)];
IO2_idx = [idx2',IO2.GRFs.right_foot(1:2),IO2.GRFs.left_foot(1:2)];

%%

Q_bounds = get_default_bounds_dummy_motion(coord_names);
Q_bounds(1,1) = -10;
Q_bounds(2,1) = 10;
Q_bounds(1,2) = 0;
Q_bounds(1,3) = 0.5;
Q_bounds(2,2) = 2;
Q_bounds(2,3) = 1;
Q_scale = diff(Q_bounds);
Qs1 = lhsdesign(1,n_coord)'.*Q_scale' + Q_bounds(1,:)';
Qs1(IO1.jointi.rotations,:) = Qs1(IO1.jointi.rotations,:)*pi/180;

Qs1(idx1_l) = Qs1(idx1_r); % symmetric motion

Qdots1 = (lhsdesign(n_coord,1)-0.5)*5*0;
Qddots1 = (lhsdesign(n_coord,1)-0.5)*10*0;

Qs2 = zeros(n_coord2,1);
Qs2(idx2,1) = Qs1(idx1);
Qdots2(idx2,1) = Qdots1(idx1);
Qddots2(idx2,1) = Qddots1(idx1);

QsQdots1 = zeros(n_coord*2,1);
QsQdots1(1:2:end) = Qs1(:);
QsQdots1(2:2:end) = Qdots1(:);

QsQdots2 = zeros(n_coord2*2,1);
QsQdots2(1:2:end) = Qs2(:);
QsQdots2(2:2:end) = Qdots2(:);

%%

res1 = F1([QsQdots1;Qddots1]);
res1 = full(res1);
res1 = res1(IO1_idx);

res2 = F2([QsQdots2;Qddots2]);
res2_0 = full(res2);
% res2 = res2_0;
res2 = res2_0(IO2_idx);

res2_reo = res2;
% res2_reo(idx1) = res2(idx2);

%% 
diff = res2_reo - res1;

max(abs(diff))

iserr = nan(size(res1));
for i=1:length(iserr)
    iserr(i) = abs(res2_reo(i) - res1(i)) > eps(res1(i));
end

find(iserr)

res_all = [res1,res2_reo,diff,diff./res1*100];

%%

% Qs_MX = MX.sym('in',n_coord,1);
% Qdots_MX = MX.sym('in',n_coord,1);
% Qddots_MX = MX.sym('in',n_coord,1);
% 
% QsQdots_MX = MX.zeros(n_coord*2,1);
% QsQdots_MX(1:2:end) = Qs_MX(:);
% QsQdots_MX(2:2:end) = Qdots_MX(:);
% 
% 
% res1_MX = F1([QsQdots_MX;Qddots_MX]);
% res2_MX = F2([QsQdots_MX;Qddots_MX]);

%%

% figure
% subplot(121)
% spy(jacobian(res1_MX,[Qs_MX;Qdots_MX;Qddots_MX]));
% xlabel('Qs, Qdots, Qddots')
% title('old external function')
% xline(n_coord+0.5,'--k')
% xline(2*n_coord+0.5,'--k')
% cnt = n_coord;
% yline(cnt+0.5,'--k')
% cnt = cnt + length(fields(IO1.origin))*3;
% yline(cnt+0.5,'--k')
% cnt = cnt + length(fields(IO1.GRFs))*3;
% yline(cnt+0.5,'--k')
% cnt = cnt + length(fields(IO1.GRMs))*3;
% yline(cnt+0.5,'--k')
% cnt = cnt + length(fields(IO1.P_contact_deformation_y));
% yline(cnt+0.5,'--k')
% ylim([0,cnt+1])
% xlim([0,n_coord*3+1])
% 
% 
% subplot(122)
% spy(jacobian(res2_MX,[Qs_MX;Qdots_MX;Qddots_MX]));
% axis tight
% xlabel('Qs, Qdots, Qddots')
% title('new external function')
% xline(n_coord+0.5,'--k')
% xline(2*n_coord+0.5,'--k')
% cnt = n_coord;
% yline(cnt+0.5,'--k')
% cnt = cnt + length(fields(IO1.origin))*3;
% yline(cnt+0.5,'--k')
% cnt = cnt + length(fields(IO1.GRFs))*3;
% yline(cnt+0.5,'--k')
% cnt = cnt + length(fields(IO1.GRMs))*3;
% yline(cnt+0.5,'--k')
% cnt = cnt + length(fields(IO1.P_contact_deformation_y));
% yline(cnt+0.5,'--k')
% ylim([0,cnt+1])
% xlim([0,n_coord*3+1])
% 
% sgtitle('Full Jacobian sparsity')


% figure
% subplot(3,2,[1,3])
% cnt = n_coord + length(fields(IO1.origin))*3;
% spy(jacobian(res1_MX(1:cnt),[Qs_MX;Qdots_MX;Qddots_MX]));
% xlabel('Qs, Qdots, Qddots')
% ylabel('positions    ID torques')
% title('old external function')
% xline(n_coord+0.5,'-','Color',[1,1,1]*0.6)
% xline(2*n_coord+0.5,'-','Color',[1,1,1]*0.6)
% cnt = n_coord;
% yline(cnt+0.5,'-','Color',[1,1,1]*0.6)
% cnt = cnt + length(fields(IO1.origin))*3;
% ylim([0,cnt+1])
% xlim([0,n_coord*3+1])
% set(gca,'XTick',[1,n_coord,n_coord*2,n_coord*3])
% set(gca,'YTick',[1,n_coord,cnt])
% 
% subplot(3,2,[2,4])
% cnt = n_coord + length(fields(IO1.origin))*3;
% spy(jacobian(res2_MX(1:cnt),[Qs_MX;Qdots_MX;Qddots_MX]));
% axis tight
% xlabel('Qs, Qdots, Qddots')
% ylabel('positions    ID torques')
% title('new external function')
% xline(n_coord+0.5,'-','Color',[1,1,1]*0.6)
% xline(2*n_coord+0.5,'-','Color',[1,1,1]*0.6)
% cnt = n_coord;
% yline(cnt+0.5,'-','Color',[1,1,1]*0.6)
% cnt = cnt + length(fields(IO1.origin))*3;
% ylim([0,cnt+1])
% xlim([0,n_coord*3+1])
% set(gca,'XTick',[n_coord,n_coord*2,n_coord*3])
% set(gca,'YTick',[n_coord,cnt])
% 
% sgtitle('OCP Jacobian sparsity')





