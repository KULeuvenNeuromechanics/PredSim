clear
close all
clc

%% input data
data_folder = 'C:\Users\u0150099\OneDrive - KU Leuven\PhD\AFO project\2023_10_03';


% ankle measurements
data_files{1} = '2023-10-03_11-17-15_2023_10_03_left_ankle.csv';
time_windows{1} = [2.5,13]; % s

data_files{2} = '2023-10-03_14-42-07_2023_10_03_left_ankle.csv';
time_windows{2} = [2.5,13]; % s

cutoff_freq = 2; % Hz, lowpass

modelfun = @(coeff, x) coeff(1) + coeff(2).*x(:,1) + coeff(3).*x(:,2);
coeff_0 = [-15,500,100];




figure
hold on
plot(q_raw,T_raw,'DisplayName','raw')
plot(q_filt,T_filt,'DisplayName','filtered')
legend('Location','best')


% MTP measurements
% data_files = '2023-10-03_14-46-13_2023_10_03_left_tip.csv';
% time_window = [2.8,11.8]; % s
% 
% cutoff_freq = 2; % Hz, lowpass
% 
% modelfun = @(coeff, x) ( coeff(1) + coeff(2)*exp((coeff(3)+x(:,1)).*coeff(4)) ).*(coeff(5).*x(:,2) + 1);
% coeff_0 = [-0.2,0.2,-0.04/28,28,2];







%%
function [q_filt, qdot_filt, T_filt, varargout] = getFilteredData(data_file, options)

arguments
    data_file char
    options.sampling_freq (1,1) double = 100; % Hz
    options.data_folder char {isfolder} = pwd;
    options.time_window (1,2) double = [0,inf]; % s
    options.cutoff_freq (1,1) double;
    options.ButterWorthFilterOrder (1,1) double {isinteger} = 6;
end

sampling_freq = options.sampling_freq;
data_folder = options.data_folder;
time_window = options.time_window;
cutoff_freq = options.cutoff_freq;
filter_order = options.ButterWorthFilterOrder;


% get matrix with data
dat = importdata(fullfile(data_folder,data_file));
if isstruct(dat)
    dat = dat.data;
end

% read data
time = dat(:,1)/sampling_freq; % [s]
q_raw = dat(:,2)*pi/180; % [rad]
T_raw = dat(:,3); % [Nm]


% only keep data inside time window
idx = find(time>time_window(1) & time<time_window(2));
time = time(idx);
q_raw = q_raw(idx);
T_raw = T_raw(idx);

% use spline to estimate velocity
qdot_raw = fnval(fnder(spline(time,q_raw)),time);

% optionally return raw data
if nargout > 3
    varargout{1} = q_raw;
end
if nargout > 4
    varargout{2} = qdot_raw;
end
if nargout > 5
    varargout{3} = T_raw;
end




q_f = fftshift(fft(q_raw));
T_f = fftshift(fft(T_raw));


figure
Tq = T_f./q_f;
f = sampling_freq/2*linspace(-1,1,length(Tq));
semilogx(f,db(abs(Tq)))
hold on
xline(cutoff_freq,'-k')


[b,a] = butter(filter_order,cutoff_freq/sampling_freq*2);

% figure
% freqz(b,a)

q_filt = filter(b,a,q_raw);
T_filt = filter(b,a,T_raw);

qdot_filt = fnval(fnder(spline(time,q_filt)),time);




end

%%
function [] = fitFunctionToData()
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


end
