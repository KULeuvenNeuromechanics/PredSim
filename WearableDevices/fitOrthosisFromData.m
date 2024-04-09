% --------------------------------------------------------------------------
% fitOrthosisFromData
%   Script to fit a function describing the orthosis based on measured
%   data. Orthosis function calculates a moment from the bending angle and the 
%   bending angular velocity.
% 
% note: this is still a draft version, but feedback is welcome.
% 
% Original author: Lars D'Hondt
% Original date: January 2024
% --------------------------------------------------------------------------



clear
close all
clc

%% input data

data_folder = 'C:\Users\u0150099\Documents\orthosistester\test 2';
time_window = {[7,33],[3,16],[2.5,7]};
cutoff_freq = [2,3,5];


% data_folder = 'C:\Users\u0150099\OneDrive - KU Leuven\PhD\AFO project\2024_01_17';
% time_window = {[2.5,16],[2.5,25],[2.5,6.6]};
% cutoff_freq = [2,2,5];

data_files = dir(fullfile(data_folder,'*_ankle.csv'));

% modelfun = @(coeff, x) coeff(1).*( x(:,1) - coeff(2) ) + coeff(3).*x(:,2);
% coeff_0 = [-15,500,100];

modelfun = @(coeff, x) (coeff(1) + coeff(2).*( x(:,1)-coeff(3) ).^2) .*...
    (1 + coeff(4).*x(:,2)) - coeff(5)*coeff(4).*x(:,2);

% modelfun = @(coeff, x) (coeff(1) + coeff(2).*( x(:,1)-coeff(3) ).^2);


f1=figure;
tiledlayout('flow')
p1 = [];
clrs = parula(1+2*length(data_files));

f2 = figure;

q_all = [];
qdot_all = [];
T_all = [];

for i=1:length(data_files)

    [q_filt, qdot_filt, T_filt, time, q, qdot, T] = getFilteredData(data_files(i).name,...
        'data_folder',data_folder, 'cutoff_freq',cutoff_freq(i), 'time_window',time_window{i});

    q_all = [q_all; q_filt];
    qdot_all = [qdot_all; qdot_filt];
    T_all = [T_all; T_filt];

    figure(f1)
    nexttile(1)
    hold on
    p1(end+1) = plot(time,q*180/pi,'Color',clrs(2*i-1,:),'DisplayName',sprintf("%i: raw",i));
    p1(end+1) = plot(time,q_filt*180/pi,'Color',clrs(2*i,:),'LineWidth',1,...
        'DisplayName',sprintf("%i: %ith order low-pass with cut-off %d Hz",i,6,cutoff_freq(i)));
    xlabel('time [s]')
    ylabel('angle [°]')

    nexttile(2)
    hold on
    plot(time,qdot,'Color',clrs(2*i-1,:))
    plot(time,qdot_filt,'Color',clrs(2*i,:),'LineWidth',1)
    xlabel('time [s]')
    ylabel('velocity [rad/s]')
    
    nexttile(3)
    hold on
    plot(time,T','Color',clrs(2*i-1,:))
    plot(time,T_filt,'Color',clrs(2*i,:),'LineWidth',1)
    xlabel('time [s]')
    ylabel('torque [Nm]')

    nexttile(4)
    hold on
    plot(q*180/pi,T,'.-','Color',clrs(2*i-1,:))
    plot(q_filt*180/pi,T_filt,'Color',clrs(2*i,:),'LineWidth',1)
    xlabel('angle [°]')
    ylabel('torque [Nm]')
    

    figure(f2)
    plot3(q_filt,qdot_filt,T_filt,'Color',clrs(2*i,:))
    hold on
    xlabel('angle [rad]')
    ylabel('velocity [rad/s]')
    zlabel('torque [Nm]')
    grid on

end
lg = legend(p1);
lg.Orientation = 'horizontal';
lg.NumColumns = 2;
lg.Layout.Tile = 'south';


%%
coeff_0 = [3,-30,0.5,-0.1,10];
% coeff_0 = [3,-30,0.5];

[f_fit, T_fit] = fitFunctionToData(modelfun, coeff_0, q_all, qdot_all, T_all, true);

figure(f1)
nexttile(4)
plot(q_all*180/pi,T_fit)

%%

figure
tiledlayout('flow')
nexttile
T_diff = T_all - T_fit;
plot(q_all*180/pi,T_diff,'.')
xlabel('q')
ylabel('T diff')

nexttile
plot(qdot_all*180/pi,T_diff,'.')
xlabel('qdot')
ylabel('T diff')

nexttile
plot(q_all*180/pi,T_diff./qdot_all,'.')
xlabel('q')
ylabel('T diff / qdot')

%%

% figure
% hold on
% plot(q_raw,T_raw,'DisplayName','raw')
% plot(q_filt,T_filt,'DisplayName','filtered')
% legend('Location','best')


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

% tmp = time(min(q_raw(idx([1:3]))) <= q_raw &...
%     q_raw <= max(q_raw(idx([1:3]))) &...
%     time<time_window(2));
% time_window(2) = tmp(end);

tmp = time(min(q_raw(idx([end-2:end]))) <= q_raw &...
    q_raw <= max(q_raw(idx([end-2:end]))) &...
    time>time_window(1));
time_window(1) = tmp(1);

idx = find(time>time_window(1) & time<time_window(2));
q_offset = mean( q_raw(time<time_window(1) & time<2.5) );
T_offset = mean( T_raw(time<time_window(1) & time<2.5) );
time = time(idx);
q_raw = q_raw(idx);
T_raw = T_raw(idx);


% use spline to estimate velocity
qdot_raw = fnval(fnder(spline(time,q_raw)),time);

% optionally return raw data
if nargout > 3
    varargout{1} = time;
end
if nargout > 4
    varargout{2} = q_raw;
end
if nargout > 5
    varargout{3} = qdot_raw;
end
if nargout > 6
    varargout{4} = T_raw;
end




q_f = fftshift(fft(q_raw));
T_f = fftshift(fft(T_raw));
Tq = T_f./q_f;

figure
f = sampling_freq/2*linspace(-1,1,length(Tq));
semilogx(f(f>0),db(abs(Tq(f>0))))
hold on
xline(cutoff_freq,'-k')


[b,a] = butter(filter_order,cutoff_freq/sampling_freq*2);

% figure
% freqz(b,a)

q_filt = filter(b,a,q_raw);
T_filt = filter(b,a,T_raw);

qdot_filt = fnval(fnder(spline(time,q_filt)),time);

q_filt_f = fftshift(fft(q_filt));
T_filt_f = fftshift(fft(T_filt));
Tq_filt = T_filt_f./q_filt_f;
semilogx(f(f>0),db(abs(Tq_filt(f>0))))


end

%%
function [f_fittedCurve, T_fit] = fitFunctionToData(modelfun, coeff_0, q_filt, qdot_filt, T_filt, manual)

if manual
    coeff_sol = coeff_0;

else
    % data points
    x_fit = [q_filt, qdot_filt];
    y_fit = T_filt;
    % fit the model
    mdl = fitnlm(x_fit,y_fit,modelfun,coeff_0);
    % get coefficient values
    coeff_sol = table2array(mdl.Coefficients(:,1));

end

% fill in the fitted coefficients
f_fittedCurve = @(q) modelfun(coeff_sol,q);

T_fit = f_fittedCurve([q_filt,qdot_filt]);


%%
% figure
% tiledlayout('flow')
% nexttile
% hold on
% plot(time,q)
% plot(time,q_filt)
% 
% nexttile
% hold on
% plot(time,qdot_ankle)
% plot(time,qdot_filt)
% 
% nexttile
% hold on
% plot(time,T)
% plot(time,T_filt)
% 
% 
% plot(time,T_fit)


end
