% This script generates a series of plots showing bounds and initial guess
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
suptitle_Fontsize = 16;
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
s = suptitle('Qs and Qdots');
set(s,'Fontsize',suptitle_Fontsize)
figure()
for i = 1:size(bounds.Qdotdots.lower,2)
    subplot(6,6,i)
    plot([1,N],[bounds.Qdotdots.upper(:,i),bounds.Qdotdots.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.Qdotdots.lower(:,i),bounds.Qdotdots.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.Qdotdots(:,i),'k','linewidth',2);
end
s = suptitle('Time derivative of Qdots');
set(s,'Fontsize',suptitle_Fontsize)
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
s = suptitle('Muscle activations');
set(s,'Fontsize',suptitle_Fontsize)
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
s = suptitle('Time derivative of muscle activations');
set(s,'Fontsize',suptitle_Fontsize)
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
s = suptitle('Muscle-tendon forces');
set(s,'Fontsize',suptitle_Fontsize)
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
s = suptitle('Time derivative of muscle-tendon forces');
set(s,'Fontsize',suptitle_Fontsize)
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
s = suptitle('Arm activations');
set(s,'Fontsize',suptitle_Fontsize)
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
s = suptitle('Arm excitations');
set(s,'Fontsize',suptitle_Fontsize)

figure()
for i = 1:nq.mtp
    subplot(3,3,i)
    plot([1,N],[bounds.a_mtp.upper(:,i),bounds.a_mtp.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.a_mtp.lower(:,i),bounds.a_mtp.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.a_mtp(:,i),'k','linewidth',2);
end
s = suptitle('Mtp activations');
set(s,'Fontsize',suptitle_Fontsize)
figure()
for i = 1:nq.mtp
    subplot(3,3,i)
    plot([1,N],[bounds.e_mtp.upper(:,i),bounds.e_mtp.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.e_mtp.lower(:,i),bounds.e_mtp.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.e_mtp(:,i),'k','linewidth',2);
end
s = suptitle('Mtp excitations');
set(s,'Fontsize',suptitle_Fontsize)
