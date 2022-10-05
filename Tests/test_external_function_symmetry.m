
% Lars D'Hondt
% 17/May/2022


clear
clc
import casadi.*

% add repository to workspace
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/PreProcessing'])


name_1 = 'subject1_2D_3D';
% name_1 = 'PredSim_2D';
F1  = external('F',fullfile(pathRepo,'Subjects',name_1,['F_' name_1 '.dll'])); 

load(fullfile(pathRepo,'Subjects',name_1,['F_' name_1 '_IO.mat'])); 
IO1 = IO;

% F1 = external('F','C:\Users\u0150099\OneDrive - KU Leuven\PhD\algodiff\predictiveSimulations_2D\ExternalFunctions\PredSim_2D_pp.dll');
% IO1.coordi.pelvis_tilt  = 1; 
% IO1.coordi.pelvis_tx    = 2;
% IO1.coordi.pelvis_ty    = 3;
% IO1.coordi.hip_flexion_l        = 4;
% IO1.coordi.hip_flexion_r        = 5;
% IO1.coordi.knee_angle_l       = 6;
% IO1.coordi.knee_angle_r       = 7;
% IO1.coordi.ankle_angle_l      = 8;
% IO1.coordi.ankle_angle_r      = 9;
% IO1.coordi.lumbar_extension    = 10;
% IO1.GRFs.right_foot = [11,12];
% IO1.GRFs.left_foot = [13,14];

%%

coord_names = fieldnames(IO1.coordi);
n_coord = length(coord_names);

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

Qdots1 = (lhsdesign(n_coord,1)-0.5)*5;
Qddots1 = (lhsdesign(n_coord,1)-0.5)*10;

Qs2(:,1) = Qs1(:);
Qdots2(:,1) = Qdots1(:);
Qddots2(:,1) = Qddots1(:);

invA = [IO1.jointi.leg_r,IO1.jointi.leg_l];
invB = [IO1.jointi.leg_l,IO1.jointi.leg_r];

Qs2(invA) = Qs1(invB);
Qdots2(invA) = Qdots1(invB);
Qddots2(invA) = Qddots1(invB);

QsQdots1 = zeros(n_coord*2,1);
QsQdots1(1:2:end) = Qs1(:);
QsQdots1(2:2:end) = Qdots1(:);

QsQdots2 = zeros(n_coord*2,1);
QsQdots2(1:2:end) = Qs2(:);
QsQdots2(2:2:end) = Qdots2(:);

%%

res1 = F1([QsQdots1;Qddots1]);
res1 = full(res1);

res2 = F1([QsQdots2;Qddots2]);
res2 = full(res2);
res2(invA) = res2(invB);
% idx_inv_grfA = [IO1.GRFs.right_foot,IO1.GRFs.left_foot,IO1.GRFs.contact_sphere_1,IO1.GRFs.contact_sphere_2,IO1.GRFs.contact_sphere_3,IO1.GRFs.contact_sphere_4];
% idx_inv_grfB = [IO1.GRFs.left_foot,IO1.GRFs.right_foot,IO1.GRFs.contact_sphere_3,IO1.GRFs.contact_sphere_4,IO1.GRFs.contact_sphere_1,IO1.GRFs.contact_sphere_2];
idx_inv_grfA = [IO1.GRFs.right_foot,IO1.GRFs.left_foot];
idx_inv_grfB = [IO1.GRFs.left_foot,IO1.GRFs.right_foot];

res2(idx_inv_grfA) = res2(idx_inv_grfB);

diff = res1 - res2;

is_symm = diff <= eps(res1);
is_symm = is_symm(1:length(coord_names));
coord_names(~is_symm);


%%


QsQdots1 = zeros(n_coord*2,1);
QsQdots3 = zeros(n_coord*2,1);

Qs3 = Qs1;
Qs3(IO.coordi.pelvis_tx) = Qs3(IO.coordi.pelvis_tx) + 0.2;

QsQdots1(1:2:end) = Qs1(:);
QsQdots1(2:2:end) = Qdots1(:);

QsQdots3(1:2:end) = Qs3(:);
QsQdots3(2:2:end) = Qdots1(:);

res1 = F1([QsQdots1;Qddots1]);
res1 = full(res1);

res3 = F1([QsQdots3;Qddots1]);
res3 = full(res3);

diff = res3 - res1;
max(abs(diff))

res_all = [res1,res3,diff,diff./res1*100];
