function [IC1i_c,IC1i_s,HS1,threshold] = getHeelstrikeSimulation(GRFk_opt,N,threshold_init)
% --------------------------------------------------------------------------
% getHeelstrikeSimulation
%  To reconstruct the full gait cycle from the simulated half gait cycle
%   
% INPUT:
%   - GRFk_opt -
%   * Ground reaction forces. Columns are x,y,z for right foot, then x,y,z
%   for left foot.
%
%   - N -
%   * Number of mesh intervals
%
%   - threshold_init -
%   * initial value of the threshold for vertical GRF (in N)
%
% OUTPUT:
%   - IC1i_c -
%   * index of initial contact for controls
%
%   - IC1i_s -
%   * index of initial contact for states
%
%   - HS1 -
%   * which foot makes inital contact 'r' or 'l'
% 
%   - threshold -
%   * value of the threshold for vertical GRF (in N)
%
% Original author: Maarten Afschrift
% Original date: 29/Nov/2021
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

% Identify heel strike
threshold = threshold_init; % there is foot-ground contact above the threshold
if exist('HS1','var')
    clear HS1
end

% increase threshold untill you have at least one frame above the threshold
nFramesBelow= sum(GRFk_opt(:,2)<threshold);
while nFramesBelow == 0
    threshold = threshold + 1;
    nFramesBelow= sum(GRFk_opt(:,2)<threshold);
end
if threshold <100
    phase_tran_tgridi = find(GRFk_opt(:,2)<threshold,1,'last');
else
    % heelstrike is in between left and right leg simulation
    if threshold >100 && GRFk_opt(end,5)<20
        phase_tran_tgridi = find(GRFk_opt(:,2)<threshold,1,'last');
    else
        % heelstrike is on the left leg
        phase_tran_tgridi =[];
        threshold = threshold_init;
    end
end
if ~isempty(phase_tran_tgridi)
    if phase_tran_tgridi == N
        temp_idx = find(GRFk_opt(:,2)>threshold,1,'first');
        if ~isempty(temp_idx)
            if temp_idx-1 ~= 0 && ...
                    find(GRFk_opt(temp_idx-1,2)<threshold)
                phase_tran_tgridi_t = temp_idx;
                IC1i = phase_tran_tgridi_t;
                HS1 = 'r';
            end
        else
            IC1i = phase_tran_tgridi + 1;
            HS1 = 'r';
        end
    else
        IC1i = phase_tran_tgridi + 1;
        HS1 = 'r';
    end
end
if ~exist('HS1','var')
    % Check if heel strike is on the left side
    phase_tran_tgridi = find(GRFk_opt(:,5)<threshold,1,'last');
    if phase_tran_tgridi == N
        temp_idx = find(GRFk_opt(:,5)>threshold,1,'first');
        if ~isempty(temp_idx)
            if temp_idx-1 ~= 0 && ...
                    find(GRFk_opt(temp_idx-1,5)<threshold)
                phase_tran_tgridi_t = temp_idx;
                IC1i = phase_tran_tgridi_t;
                HS1 = 'l';
            else
                IC1i = phase_tran_tgridi + 1;
                HS1 = 'l';
            end
        else
            IC1i = phase_tran_tgridi + 1;
            HS1 = 'l';
        end
    else
        IC1i = phase_tran_tgridi + 1;
        HS1 = 'l';
    end
end

% GRFk_opt is at mesh points starting from k=2, we thus add 1 to IC1i
% for the states
if phase_tran_tgridi ~= N
    IC1i_c = IC1i;
    IC1i_s = IC1i + 1;
else
    IC1i_c = IC1i;
    IC1i_s = IC1i;
end


end

