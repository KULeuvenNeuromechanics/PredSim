% This function returns the Qs (radian) given the path to the inverse
% kinematics file and the joints of interest.
%
% Author: Antoine Falisse
% Date: 12/19/2018
% Adapted: Lars D'Hondt
% Date: 30 nov 2021
%--------------------------------------------------------------------------

function Qs_IK = getIK(path_IK,model_info)

% Get the names of the coordinates
coordinate_names = fieldnames(model_info.ExtFunIO.coordi);

% Load the reference date
Qsall = importdata(path_IK);

%% Build the struct with inverse kinematics
% First column has the timestamps
Qs_IK.time = Qsall.data(:,strcmp(Qsall.colheaders,{'time'}));
Qs_IK.all(:,1) = Qs_IK.time;
Qs_IK.colheaders{1,1} = 'time';

% Add the IK of every joint, in the defined order
for i = 1:length(coordinate_names)
    coordinate = coordinate_names{i};
    idx_coord_Qsall = find(strcmp(Qsall.colheaders,coordinate),1,"first");
    if isempty(idx_coord_Qsall)
        Qs_IK.(coordinate) = zeros(size(Qs_IK.time));
    else
        if model_info.jointi.translation(:) == model_info.coordi.(coordinate)
            Qs_IK.(coordinate) = Qsall.data(:,strcmp(Qsall.colheaders,coordinate));
        else
            Qs_IK.(coordinate) = Qsall.data(:,strcmp(Qsall.colheaders,coordinate)).*(pi/180);
        end
    end
    Qs_IK.all(:,model_info.coordi.(coordinate)+1) = Qs_IK.(coordinate);
    Qs_IK.colheaders(1,model_info.coordi.(coordinate)+1) = coordinate;
end

% Low-pass filter
order = 4;
cutoff_low = 6;
fs=1/mean(diff(Qs_IK.all(:,1)));
[af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
Qs_IK.allfilt = Qs_IK.all;
Qs_IK.allfilt(:,2:end) = filtfilt(af,bf,Qs_IK.allfilt(:,2:end));

end
