clear
clc

import casadi.*

N = 200;
q = linspace(-100,10,N);

f_lMT = Function.load('Vitruvian_Man_f_lMT_vMT_dM_poly_3_9');

load('F_Vitruvian_Man_IO.mat')


lMT = nan(N,5);
dM = nan(N,5);

for i=1:N
    qi = zeros(12,1);
    qi(IO.coordi.knee_angle_r,1) = q(i)*pi/180;
    
    
    [lMTi,~,dMi] = f_lMT(qi,qi*0);
    dMi = full(dMi(:,IO.coordi.knee_angle_r));
    idx = find(dMi~=0);
    
    lMT(i,:) = full(lMTi(idx));
    dM(i,:) = dMi(idx);
    
end

% mus_names = model_info.muscle_info.muscle_names(idx);

mus_names = {'hamstrings_r','bifemsh_r','rect_fem_r','vasti_r','gastroc_r'};


%%

figure
tiledlayout('flow')
for i=1:5
    nexttile(i)
    hold on
    plot(q,dM(:,i))
    ylabel('dM (m)')
    xlabel('q knee (°)')
    xline(-90,'Color','k')
    xline(0,'Color','k')
    title(replace(mus_names{i},'_',' '))
    
    
end
