function Qs_IK = getIK(path_IK,model_info)
% --------------------------------------------------------------------------
% getIK
%   This function returns the coordinate positions (radian) given the path 
%   to the inverse kinematics file and the joints of interest.
%   
% INPUT:
%   - path_IK -
%   * path to the inverse kinematics file
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - Qs_IK -
%   * struct with inverse kinematics results (time, coordinates and labels)
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
% --------------------------------------------------------------------------
% This file is part of PredSim.
% 
% PredSim: A Framework for Rapid Predictive Simulations of Locomotion
% Copyright (c) 2026 KU Leuven
% 
% PredSim is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Affero General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% 
% PredSim is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public 
% License for more details.
% 
% You should have received a copy of the GNU Affero General Public License 
% along with PredSim. If not, see <https://www.gnu.org/licenses/>.
% --------------------------------------------------------------------------

% Get the names of the coordinates
coordinate_names = fieldnames(model_info.ExtFunIO.coordi);

% Load the reference date
Qsall = read_motionFile_v40(path_IK);

%% Build the struct with inverse kinematics
% First column has the timestamps
Qs_IK.time = Qsall.data(:,strcmp(Qsall.labels,{'time'}));
Qs_IK.all(:,1) = Qs_IK.time;
Qs_IK.colheaders{1,1} = 'time';

% Add the IK of every joint, in the defined order
for i = 1:length(coordinate_names)
    coordinate = coordinate_names{i};
    idx_coord_Qsall = find(strcmp(Qsall.labels,coordinate),1,'first');
    if isempty(idx_coord_Qsall)
        Qs_IK.(coordinate) = zeros(size(Qs_IK.time));
    else
        if strcmp(Qsall.inDeg,'yes')
            if sum(model_info.ExtFunIO.jointi.translations(:) == model_info.ExtFunIO.coordi.(coordinate))
                Qs_IK.(coordinate) = Qsall.data(:,idx_coord_Qsall);
            else
                Qs_IK.(coordinate) = Qsall.data(:,idx_coord_Qsall).*(pi/180);
            end
        else
            Qs_IK.(coordinate) = Qsall.data(:,idx_coord_Qsall);
        end

    end
    Qs_IK.all(:,model_info.ExtFunIO.coordi.(coordinate)+1) = Qs_IK.(coordinate);
    Qs_IK.colheaders{1,model_info.ExtFunIO.coordi.(coordinate)+1} = coordinate;
end

% Low-pass filter
order = 4;
cutoff_low = 6;
fs=1/mean(diff(Qs_IK.all(:,1)));
[af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
Qs_IK.allfilt = Qs_IK.all;
Qs_IK.allfilt(:,2:end) = filtfilt(af,bf,Qs_IK.allfilt(:,2:end));

end
