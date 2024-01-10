clear
close all
clc


data_folder = 'C:\Users\u0150099\OneDrive - KU Leuven\PhD\AFO project\2023_10_03';

% data_file = '2023-10-03_11-17-15_2023_10_03_left_ankle.csv';
% cutoff_freq = 2; % Hz, lowpass
% time_window = [2.5,13]; % s
% modelfun = @(coeff, x) coeff(1) + coeff(2).*x(:,1) + coeff(3).*x(:,2);
% coeff_0 = [-15,500,100];

data_file = '2023-10-03_14-46-13_2023_10_03_left_tip.csv';
cutoff_freq = 2; % Hz, lowpass
time_window = [2.8,11.8]; % s
modelfun = @(coeff, x) ( coeff(1) + coeff(2)*exp((coeff(3)+x(:,1)).*coeff(4)) ).*(coeff(5).*x(:,2) + 1);
coeff_0 = [-0.2,0.2,-0.04/28,28,2];






sampling_freq = 100; % Hz



dat = importdata(fullfile(data_folder,data_file));

if isstruct(dat)
    dat = dat.data;
end

time = dat(:,1)/sampling_freq;
q = dat(:,2)*pi/180;
T = dat(:,3);

idx = find(time>time_window(1) & time<time_window(2));
time = time(idx);
q = q(idx);
T = T(idx);

qdot_ankle = fnval(fnder(spline(time,q)),time);


q_f = fftshift(fft(q));
T_f = fftshift(fft(T));

q_f = T_f./q_f;

f = sampling_freq/2*linspace(-1,1,length(q_f));



figure
semilogx(f,db(abs(q_f)))
hold on
xline(cutoff_freq,'-k')


[b,a] = butter(6,cutoff_freq/sampling_freq*2);

% figure
% freqz(b,a)

q_filt = filter(b,a,q);
qdot_filt = fnval(fnder(spline(time,q_filt)),time);
T_filt = filter(b,a,T);






figure
hold on
plot(q,T,'DisplayName','raw')
plot(q_filt,T_filt,'DisplayName','filtered')
legend('Location','best')

%%
coeff_sol = coeff_0;

% data points
x_fit = [q_filt, qdot_filt];
y_fit = T_filt;
% fit the model
mdl = fitnlm(x_fit,y_fit,modelfun,coeff_0);
% get coefficient values
coeff_sol = table2array(mdl.Coefficients(:,1));



% fill in the fitted coefficients
f_fittedCurve = @(q) modelfun(coeff_sol,q);

T_fit = f_fittedCurve([q_filt,qdot_filt]);

plot(q_filt,T_fit,'DisplayName','fitted')



%%
figure
tiledlayout('flow')
nexttile
hold on
plot(time,q)
plot(time,q_filt)

nexttile
hold on
plot(time,qdot_ankle)
plot(time,qdot_filt)

nexttile
hold on
plot(time,T)
plot(time,T_filt)


plot(time,T_fit)

