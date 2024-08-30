
clear
close all
clc

[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '\CasadiFunctions']);
addpath([pathRepo '\PreProcessing']);



lTtilde = linspace(1,1.1,100);
lMtilde = linspace(0.5,1.6,100);
f_get_forces = @(s1,s2,s3,s4) f_get_forces_l(lTtilde,lMtilde,s1,s2,s3,s4);

% settings
tendon_stiff_vec = [0.7,1,1.2];
muscle_strength_vec = [0.7,1,1.2];
muscle_pass_stiff_shift_vec = [1.1,1,0.9];
muscle_pass_stiff_scale_vec = [0.7,1,1.2];


[~,F_muscle_a_0,F_muscle_p_0] = f_get_forces(1,1,1,1);

ylm = [0,max([1.5, muscle_strength_vec+0.2])];

f1=figure('Position',[400 500 1400 420]);
sgtitle('Effects of scale factors on normalized muscle-tendon properties')

%% plot tendon stiffness

subplot(1,4,1)
hold on
for i=1:length(tendon_stiff_vec)
    % calculate force
    [F_tendon,~,~] = f_get_forces(tendon_stiff_vec(i),1,1,1);
    % add plot
    plot(lTtilde,F_tendon,'DisplayName',['tendon_stiff = ' num2str(tendon_stiff_vec(i))])
end

xlabel('Tendon length (1/lTs)')
ylabel('Tendon force (1/FMo)')
title('Scaling tendon stiffness')
legend('Location','southoutside','Interpreter','none')
xlim([min(lTtilde),max(lTtilde)])
grid on

%% plot muscle strength

subplot(1,4,2)
hold on
plh0 = {};
lg1 = [];

for i=1:length(muscle_strength_vec)
    % calculate force
    [~,F_muscle_a,~] = f_get_forces(1,muscle_strength_vec(i),1,1);
    % add plot
    p0=plot(lMtilde,F_muscle_a,'--');
    plh0{end+1} = p0;
end
plot(lMtilde,F_muscle_p_0,'-.k')
for i=1:length(muscle_strength_vec)
    % calculate force
    [~,F_muscle_a,~] = f_get_forces(1,muscle_strength_vec(i),1,1);
    % add plot
    p1=plot(lMtilde,F_muscle_a+F_muscle_p_0,'-','DisplayName',...
        ['muscle_strength = ' num2str(muscle_strength_vec(i))],'Color',plh0{i}.Color);
    lg1(end+1) = p1;
end

xlabel('Fiber length (1/lMo)')
ylabel('Isometric fiber force (1/FMo)')
title('Scaling muscle strength')
legend(lg1,'Location','southoutside','Interpreter','none')
xlim([min(lMtilde),max(lMtilde)])
ylim(ylm);
set(gca,'YTick',sort([muscle_strength_vec,ylm,0.5]))
grid on


%% plot muscle passive stiffness shift

subplot(1,4,3)
hold on
plh0 = {};
lg1 = [];

for i=1:length(muscle_pass_stiff_scale_vec)
    % calculate force
    [~,~,F_muscle_p] = f_get_forces(1,1,1,muscle_pass_stiff_scale_vec(i));
    % add plot
    p0=plot(lMtilde,F_muscle_p,'-.');
    plh0{end+1} = p0;
end
plot(lMtilde,F_muscle_a_0,'--k')
for i=1:length(muscle_pass_stiff_shift_vec)
    % calculate force
    [~,~,F_muscle_p] = f_get_forces(1,1,1,muscle_pass_stiff_scale_vec(i));
    % add plot
    p1=plot(lMtilde,F_muscle_a_0+F_muscle_p,'-','DisplayName',...
        ['muscle_pass_stiff_scale = ' num2str(muscle_pass_stiff_scale_vec(i))],'Color',plh0{i}.Color);
    lg1(end+1) = p1;
end

xlabel('Fiber length (1/lMo)')
ylabel('Isometric fiber force (1/FMo)')
title('Scaling passive muscle stiffness')
legend(lg1,'Location','southoutside','Interpreter','none')
xlim([min(lMtilde),max(lMtilde)])
ylim(ylm);
grid on

%% plot muscle passive stiffness shift

subplot(1,4,4)
hold on
plh0 = {};
lg1 = [];

for i=1:length(muscle_pass_stiff_shift_vec)
    % calculate force
    [~,~,F_muscle_p] = f_get_forces(1,1,muscle_pass_stiff_shift_vec(i),1);
    % add plot
    p0=plot(lMtilde,F_muscle_p,'-.');
    plh0{end+1} = p0;
end
plot(lMtilde,F_muscle_a_0,'--k')
for i=1:length(muscle_pass_stiff_shift_vec)
    % calculate force
    [~,~,F_muscle_p] = f_get_forces(1,1,muscle_pass_stiff_shift_vec(i),1);
    % add plot
    p1=plot(lMtilde,F_muscle_a_0+F_muscle_p,'-','DisplayName',...
        ['muscle_pass_stiff_shift = ' num2str(muscle_pass_stiff_shift_vec(i))],'Color',plh0{i}.Color);
    lg1(end+1) = p1;
end

xlabel('Fiber length (1/lMo)')
ylabel('Isometric fiber force (1/FMo)')
title('Shifting passive muscle stiffness')
legend(lg1,'Location','southoutside','Interpreter','none')
xlim([min(lMtilde),max(lMtilde)])
ylim(ylm);
grid on


%%
exportgraphics(f1,[pwd '/fig_muscle_tendon_properties_scaling.png']);

%%
function [Fl_T, Fla_M, Flp_M] = f_get_forces_l(lTtilde,lMtilde, tendon_stiff,...
    muscle_strength,muscle_pass_stiff_shift,muscle_pass_stiff_scale)

load('Fvparam.mat','Fvparam');
load('Fpparam.mat','Fpparam');
load('Faparam.mat','Faparam');

% lTtilde = linspace(1,1.1,100);
% lMtilde = linspace(0.5,1.6,100);

[Fl_T, Fla_M, Flp_M] = f_muscle_mechanics(Fvparam,Fpparam, Faparam, lTtilde,lMtilde,...
    tendon_stiff,muscle_strength,muscle_pass_stiff_shift,muscle_pass_stiff_scale);

end

%%
function [fse, Fiso, Fpetilde] = f_muscle_mechanics(Fvparamf,Fpparamf, Faparamf,...
    lTtildef,lMtildef,tendon_stiff,strength,stiffness_shift,stiffness_scale)

tendon_stiff = 35*tendon_stiff;

% tendon force-length characteristic
fse = exp((lTtildef - 0.995).*tendon_stiff)./5 - 0.25 + getShift(tendon_stiff);

% Active muscle force-length characteristic
b11 = Faparamf(1);
b21 = Faparamf(2);
b31 = Faparamf(3);
b41 = Faparamf(4);
b12 = Faparamf(5);
b22 = Faparamf(6);
b32 = Faparamf(7);
b42 = Faparamf(8);
b13 = 0.1;
b23 = 1;
b33 = 0.5*sqrt(0.5);
b43 = 0;
num3 = lMtildef-b23;
den3 = b33+b43*lMtildef;
FMtilde3 = b13*exp(-0.5*num3.^2./den3.^2);
num1 = lMtildef-b21;
den1 = b31+b41*lMtildef;
FMtilde1 = b11*exp(-0.5*num1.^2./den1.^2);
num2 = lMtildef-b22;
den2 = b32+b42*lMtildef;
FMtilde2 = b12*exp(-0.5*num2.^2./den2.^2);
FMltilde = FMtilde1+FMtilde2+FMtilde3;
Fiso = strength.*FMltilde;


% Passive muscle force-length characteristic
e0 = 0.6;
kpe = 4;
t5 = exp(kpe * (lMtildef - stiffness_shift) / (e0/stiffness_scale));
% Passive muscle force
Fpetilde = ((t5 - 0.10e1) - Fpparamf(1)) / Fpparamf(2);


end