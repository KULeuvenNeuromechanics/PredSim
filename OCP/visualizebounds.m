% This script generates a series of plots showing bounds and initial guess
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
title_Fontsize = 16;
figure()
for i = 1:size(bounds.Qs.lower,2)
    subplot(5,7,i)
    plot([1,N],[bounds.Qs.upper(:,i),bounds.Qs.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.Qs.lower(:,i),bounds.Qs.lower(:,i)],...
        'r--','linewidth',2);
    hold on
    plot(guess.Qs(:,i),'k','linewidth',2);
%     hold on
%     plot(Qs_spline.data(:,i+1)/scaling.Qs(i),'g--','linewidth',1);
%     hold on
%     plot(q_in.data(:,i+1)/scaling.Qs(i),'k--','linewidth',1);
%     title(coordinate_names{i})
end
s = title('Qs');
set(s,'Fontsize',title_Fontsize)
figure()
for i = 1:size(bounds.QsQdots.lower,2)
    subplot(10,7,i)
    plot([1,N],[bounds.QsQdots.upper(:,i),bounds.QsQdots.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.QsQdots.lower(:,i),bounds.QsQdots.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.QsQdots(:,i),'k','linewidth',2);
end
s = title('Qs and Qdots');
set(s,'Fontsize',title_Fontsize)
figure()
for i = 1:size(bounds.Qdotdots.lower,2)
    subplot(5,7,i)
    plot([1,N],[bounds.Qdotdots.upper(:,i),bounds.Qdotdots.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.Qdotdots.lower(:,i),bounds.Qdotdots.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.Qdotdots(:,i),'k','linewidth',2);
end
s = title('Time derivative of Qdots');
set(s,'Fontsize',title_Fontsize)
figure()
for i = 1:NMuscle
    subplot(10,10,i)
    plot([1,N],[bounds.a.upper(:,i),bounds.a.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.a.lower(:,i),bounds.a.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.a(:,i),'k','linewidth',2);
end
s = title('Muscle activations');
set(s,'Fontsize',title_Fontsize)
figure()
for i = 1:NMuscle
    subplot(10,10,i)
    plot([1,N],[bounds.vA.upper(:,i),bounds.vA.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.vA.lower(:,i),bounds.vA.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.vA(:,i),'k','linewidth',2);
end
s = title('Time derivative of muscle activations');
set(s,'Fontsize',title_Fontsize)
figure()
for i = 1:NMuscle
    subplot(10,10,i)
    plot([1,N],[bounds.FTtilde.upper(:,i),bounds.FTtilde.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.FTtilde.lower(:,i),bounds.FTtilde.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.FTtilde(:,i),'k','linewidth',2);
end
s = title('Muscle-tendon forces');
set(s,'Fontsize',title_Fontsize)
figure()
for i = 1:NMuscle
    subplot(10,10,i)
    plot([1,N],[bounds.dFTtilde.upper(:,i),bounds.dFTtilde.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.dFTtilde.lower(:,i),bounds.dFTtilde.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.dFTtilde(:,i),'k','linewidth',2);
end
s = title('Time derivative of muscle-tendon forces');
set(s,'Fontsize',title_Fontsize)
figure()
for i = 1:nq.arms
    subplot(3,3,i)
    plot([1,N],[bounds.a_a.upper(:,i),bounds.a_a.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.a_a.lower(:,i),bounds.a_a.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.a_a(:,i),'k','linewidth',2);
end
s = title('Arm activations');
set(s,'Fontsize',title_Fontsize)
figure()
for i = 1:nq.arms
    subplot(3,3,i)
    plot([1,N],[bounds.e_a.upper(:,i),bounds.e_a.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.e_a.lower(:,i),bounds.e_a.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.e_a(:,i),'k','linewidth',2);
end
s = title('Arm excitations');
set(s,'Fontsize',title_Fontsize)
