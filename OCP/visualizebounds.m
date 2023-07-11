% --------------------------------------------------------------------------
% visualizebounds
%    This script generates a series of plots showing bounds and initial guess
%   
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

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
    title(model_info.ExtFunIO.coord_names.all{i})
end
s = sgtitle('Qs');
set(s,'Fontsize',title_Fontsize)
figure()
for i = 1:size(bounds.Qdots.lower,2)
    subplot(10,7,i)
    plot([1,N],[bounds.Qdots.upper(:,i),bounds.Qdots.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.Qdots.lower(:,i),bounds.Qdots.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.Qdots(:,i),'k','linewidth',2);
end
s = sgtitle('Qdots');
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
s = sgtitle('Time derivative of Qdots');
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
s = sgtitle('Muscle activations');
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
s = sgtitle('Time derivative of muscle activations');
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
s = sgtitle('Muscle-tendon forces');
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
s = sgtitle('Time derivative of muscle-tendon forces');
set(s,'Fontsize',title_Fontsize)
figure()
for i = 1:nq.torqAct
    subplot(3,3,i)
    plot([1,N],[bounds.a_a.upper(:,i),bounds.a_a.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.a_a.lower(:,i),bounds.a_a.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.a_a(:,i),'k','linewidth',2);
end
s = sgtitle('Torque actuator activations');
set(s,'Fontsize',title_Fontsize)
figure()
for i = 1:nq.torqAct
    subplot(3,3,i)
    plot([1,N],[bounds.e_a.upper(:,i),bounds.e_a.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.e_a.lower(:,i),bounds.e_a.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.e_a(:,i),'k','linewidth',2);
end
s = sgtitle('Torque actuator excitations');
set(s,'Fontsize',title_Fontsize)
