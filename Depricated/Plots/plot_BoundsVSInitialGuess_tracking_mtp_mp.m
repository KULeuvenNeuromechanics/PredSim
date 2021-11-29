% This script generates a series of plots showing bounds and initial guess
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
suptitle_Fontsize = 16;
figure()
for i = 1:size(bounds(p).p.QsQdots.lower,2)
    subplot(10,7,i)
    plot([1,N],[bounds(p).p.QsQdots.upper(:,i),bounds(p).p.QsQdots.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds(p).p.QsQdots.lower(:,i),bounds(p).p.QsQdots.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess(p).p.QsQdots(:,i),'k','linewidth',2);
end
s = suptitle('Qs and Qdots');
set(s,'Fontsize',suptitle_Fontsize)
figure()
for i = 1:size(bounds(p).p.Qdotdots.lower,2)
    subplot(6,6,i)
    plot([1,N],[bounds(p).p.Qdotdots.upper(:,i),bounds(p).p.Qdotdots.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds(p).p.Qdotdots.lower(:,i),bounds(p).p.Qdotdots.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess(p).p.Qdotdots(:,i),'k','linewidth',2);
end
s = suptitle('Time derivative of Qdots');
set(s,'Fontsize',suptitle_Fontsize)
figure()
for i = 1:NMuscle
    subplot(10,10,i)
    plot([1,N],[bounds(p).p.a.upper(:,i),bounds(p).p.a.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds(p).p.a.lower(:,i),bounds(p).p.a.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess(p).p.a(:,i),'k','linewidth',2);
end
s = suptitle('Muscle activations');
set(s,'Fontsize',suptitle_Fontsize)
figure()
for i = 1:NMuscle
    subplot(10,10,i)
    plot([1,N],[bounds(p).p.vA.upper(:,i),bounds(p).p.vA.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds(p).p.vA.lower(:,i),bounds(p).p.vA.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess(p).p.vA(:,i),'k','linewidth',2);
end
s = suptitle('Time derivative of muscle activations');
set(s,'Fontsize',suptitle_Fontsize)
figure()
for i = 1:NMuscle
    subplot(10,10,i)
    plot([1,N],[bounds(p).p.FTtilde.upper(:,i),bounds(p).p.FTtilde.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds(p).p.FTtilde.lower(:,i),bounds(p).p.FTtilde.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess(p).p.FTtilde(:,i),'k','linewidth',2);
end
s = suptitle('Muscle-tendon forces');
set(s,'Fontsize',suptitle_Fontsize)
figure()
for i = 1:NMuscle
    subplot(10,10,i)
    plot([1,N],[bounds(p).p.dFTtilde.upper(:,i),bounds(p).p.dFTtilde.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds(p).p.dFTtilde.lower(:,i),bounds(p).p.dFTtilde.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess(p).p.dFTtilde(:,i),'k','linewidth',2);
end
s = suptitle('Time derivative of muscle-tendon forces');
set(s,'Fontsize',suptitle_Fontsize)
figure()
for i = 1:nq.arms
    subplot(3,3,i)
    plot([1,N],[bounds(p).p.a_a.upper(:,i),bounds(p).p.a_a.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds(p).p.a_a.lower(:,i),bounds(p).p.a_a.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess(p).p.a_a(:,i),'k','linewidth',2);
end
s = suptitle('Arm activations');
set(s,'Fontsize',suptitle_Fontsize)
figure()
for i = 1:nq.arms
    subplot(3,3,i)
    plot([1,N],[bounds(p).p.e_a.upper(:,i),bounds(p).p.e_a.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds(p).p.e_a.lower(:,i),bounds(p).p.e_a.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess(p).p.e_a(:,i),'k','linewidth',2);
end
s = suptitle('Arm excitations');
set(s,'Fontsize',suptitle_Fontsize)
