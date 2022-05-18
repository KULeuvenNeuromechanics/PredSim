
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
name_1 = 'Fal_s1_mtp';
f1  = Function.load(fullfile(pathRepo,'Subjects',name_1,'f_lMT_vMT_dM_poly_3_9')); 
f1b  = Function.load(fullfile(pathRepo,'Subjects',name_1,'f_lMT_vMT_dM_poly_3_9_b'));
f1c  = Function.load(fullfile(pathRepo,'Subjects',name_1,'f_lMT_vMT_dM_poly_3_9_c'));

load(fullfile(pathRepo,'Subjects',name_1,['F_' name_1 '_IO.mat'])); 
IO1 = IO;

% load 2nd
name_2 = 'Fal_s1_mtp_v2';
f2  = Function.load(fullfile(pathRepo,'Subjects',name_2,'f_lMT_vMT_dM_poly_3_9')); 

load(fullfile(pathRepo,'Subjects',name_2,['F_' name_2 '_IO.mat'])); 
IO2 = IO;


%%

coord_names = fieldnames(IO1.coordi);
n_coord = length(coord_names);

idx2 = zeros(n_coord,1);
idx1 = (1:n_coord)';
for i=1:n_coord
    idx2(i) = IO2.coordi.(coord_names{i});
end


%%
Q_bounds = get_default_bounds_dummy_motion(coord_names);
Q_scale = diff(Q_bounds);
Qs1 = lhsdesign(1,n_coord)'.*Q_scale' + Q_bounds(1,:)';
Qs1(IO1.jointi.rotations,:) = Qs1(IO1.jointi.rotations,:)*pi/180;

Qdots1 = (lhsdesign(n_coord,1)-0.5)*5;

Qs2(idx2,1) = Qs1(idx1);
Qdots2(idx2,1) = Qdots1(idx1);

%%

[lMT1,vMT1,dM1] = f1(Qs1,Qdots1);
lMT1 = full(lMT1);
vMT1 = full(vMT1);
dM1 = full(dM1);

[lMT1b,vMT1b,dM1b] = f1b(Qs1,Qdots1);
lMT1b = full(lMT1b);
vMT1b = full(vMT1b);
dM1b = full(dM1b);

[lMT1c,vMT1c,dM1c] = f1c(Qs1,Qdots1);
lMT1c = full(lMT1c);
vMT1c = full(vMT1c);
dM1c = full(dM1c);

[lMT2,vMT2,dM2] = f2(Qs2,Qdots2);
lMT2 = full(lMT2);
vMT2 = full(vMT2);
dM2 = full(dM2);

dM2_reo = dM2;
dM2_reo(:,idx1) = dM2(:,idx2);

%%
% diff_lMT = lMT2 - lMT1;
% diff_vMT = vMT2 - vMT1;
% diff_dM = dM2_reo - dM1;
% 
% max(abs(diff_lMT))
% max(abs(diff_vMT))
% max(abs(diff_dM),[],'all')

%
% diff_lMT = lMT1b - lMT1;
% diff_vMT = vMT1b - vMT1;
% diff_dM = dM1b - dM1;
% 
% max(abs(diff_lMT))
% max(abs(diff_vMT))
% max(abs(diff_dM),[],'all')

%
diff_lMT = lMT1c - lMT1;
diff_vMT = vMT1c - vMT1;
diff_dM = dM1c - dM1;
diff_dM_rel = diff_dM./dM1;

diff_dM_rel2 = diff_dM_rel;
diff_dM_rel2(abs(dM1(:,:))<1e-2)=0;

max(abs(diff_lMT))
max(abs(diff_vMT))
max(abs(diff_dM),[],'all')
max(abs(diff_dM_rel),[],'all')
max(abs(diff_dM_rel2),[],'all')


%%
% figure
% subplot(131)
% spy(dM1)
% title('Reference')
% axis tight
% ylabel('muscles')
% xlabel('dofs')
% 
% subplot(132)
% spy(dM2)
% title('With dofs in default order')
% axis tight
% ylabel('muscles')
% xlabel('dofs')
% 
% subplot(133)
% spy(dM2_reo)
% title('With dofs reordered')
% axis tight
% ylabel('muscles')
% xlabel('dofs')
% 
% sgtitle('Momentarms')

