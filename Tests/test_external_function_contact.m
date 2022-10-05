

clear
clc

import casadi.*

% add repository to workspace
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/PreProcessing'])


name_1 = 'subject1_2D';
% name_1 = 'PredSim_2D';
F1  = external('F',fullfile(pathRepo,'Subjects',name_1,['F_' name_1 '.dll'])); 

load(fullfile(pathRepo,'Subjects',name_1,['F_' name_1 '_IO.mat'])); 
IO1 = IO;


coord_names = fieldnames(IO1.coordi);
n_coord = length(coord_names);

Qs = zeros(n_coord,1);
Qdots = Qs;
Qddots = Qs;

Qs(IO.coordi.pelvis_ty) = 0.95;

T1 = calcID(Qs,Qdots,Qddots,F1);

x = linspace(-10,30,100);

Tall = zeros(length(T1),length(x));

for i=1:length(x)
    Qs(IO.coordi.pelvis_tx) = x(i);
    T2 = calcID(Qs,Qdots,Qddots,F1);
    Tall(:,i) = T2;
end

% diff2 = T2-T1;
% disp(diff2(1:n_coord))

figure
tiledlayout('flow')

gcf;

nexttile(2)
hold on
plot(x,Tall(IO1.GRFs.right_foot(1),:),'.')
ylabel('right GRF x')
xlabel('pelvis tx (m)')
title('right')
nexttile(4)
hold on
plot(x,Tall(IO1.GRFs.right_foot(2),:),'.')
ylabel('right GRF y')
xlabel('pelvis tx (m)')
nexttile(6)
hold on
plot(x,Tall(IO1.GRFs.right_foot(3),:),'.')
ylabel('right GRF z')
xlabel('pelvis tx (m)')

nexttile(1)
hold on
plot(x,Tall(IO1.GRFs.left_foot(1),:),'.')
ylabel('left GRF x')
xlabel('pelvis tx (m)')
title('left')
nexttile(3)
hold on
plot(x,Tall(IO1.GRFs.left_foot(2),:),'.')
ylabel('left GRF y')
xlabel('pelvis tx (m)')
nexttile(5)
hold on
plot(x,Tall(IO1.GRFs.left_foot(3),:),'.')
ylabel('left GRF z')
xlabel('pelvis tx (m)')

figure
for i=1:4
    idxsi = IO1.GRFs.(['contact_sphere_' num2str(i)]);
    for j=1:3
        subplot(3,4,i+4*(j-1))
        hold on
        plot(x,Tall(idxsi(j),:),'.')
        if j==1
            title(['contact_sphere_' num2str(i)],'Interpreter','none')
        end
    end
    xlabel('pelvis tx (m)')
end
subplot(3,4,1)
ylabel('GRF x')
subplot(3,4,5)
ylabel('GRF y')
subplot(3,4,9)
ylabel('GRF z')

%%
function [Ts] = calcID(qs,dqs,ddqs,Fext)
import casadi.*
qsdqs = zeros(length(qs)*2,1);
qsdqs(1:2:end) = qs;
qsdqs(2:2:end) = dqs;

res = Fext([qsdqs; ddqs]);
Ts = full(res);

end