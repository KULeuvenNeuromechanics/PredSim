% This function returns the Qs (radian) given the path to the inverse
% kinematics file and the joints of interest.
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function Qs = getIK_original(pathIK,joints)

Qsall = read_motionFile_v40(pathIK);

Qs.time = Qsall.data(:,strcmp(Qsall.labels,{'time'}));
Qs.all(:,1) = Qs.time;
Qs.colheaders{1,1} = 'time';
count = 1;
for i = 1:length(joints)
    count = count + 1;
    if strcmp(joints{i},'pelvis_tx') || strcmp(joints{i},'pelvis_ty') || strcmp(joints{i},'pelvis_tz')
        Qs.(joints{i}) = Qsall.data(:,strcmp(Qsall.labels,joints{i}));
    else
        Qs.(joints{i}) = Qsall.data(:,strcmp(Qsall.labels,joints{i})).*(pi/180);
    end
    Qs.all(:,count) = Qs.(joints{i});
    Qs.colheaders(1,count) = Qsall.labels(1,strcmp(Qsall.labels,joints{i}));
end
% Low-pass filter
order = 4;
cutoff_low = 6;
fs=1/mean(diff(Qs.all(:,1)));
[af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
Qs.allfilt = Qs.all;
Qs.allfilt(:,2:end) = filtfilt(af,bf,Qs.allfilt(:,2:end));

end
