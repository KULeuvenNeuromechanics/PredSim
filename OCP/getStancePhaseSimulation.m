function [idx_GC,idx_init_contact,HS1,threshold] = getStancePhaseSimulation(GRFk_opt,N,threshold_init)
% --------------------------------------------------------------------------
% getStancePhaseSimulation
%   To reconstruct the full gait cycle, starting at right heel strike, from 
%   the simulated (half) gait cycle
%   
% INPUT:
%   - GRFk_opt -
%   * Ground reaction forces. Columns are x,y,z for right foot, then x,y,z
%   for left foot. (in N)
%
%   - N -
%   * Number of mesh intervals
%
%   - threshold_init -
%   * initial value of the threshold for vertical GRF (in N)
%
% OUTPUT:
%   - idx_GC -
%   * indices of mesh points such that 1st point is initial contact
%   GRFk_opt_GC = GRFk_opt(idx_GC,:);
%
%   - HS1 -
%   * which foot makes inital contact 'r' or 'l'
% 
%   - threshold -
%   * value of the threshold for vertical GRF (in N)
%
% Original author: Lars D'Hondt
% Original date: 7/March/2023
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

% only use vertical
GRFy = GRFk_opt(:,[2,5]);

%% increase threshold untill you have at least one frame above the threshold
threshold = threshold_init; 
nFramesBelow= sum(GRFy < threshold,"all");
while nFramesBelow == 0
    threshold = threshold + 1;
    nFramesBelow= sum(GRFy < threshold,"all");
    if threshold-threshold_init > 200
        error('Vertical ground reaction forces are all way above heelstrike threshold.') % just in case
    end
end

%% We define start of gait cycle as initial contact of right foot. 
% In the (unlikely) case that there is no ground contact on the right foot, 
% we consider the left foot initial contact as starting point.
is_stance = GRFy(:,1) >= threshold;
if sum(is_stance) == 0
    is_stance = GRFy(:,2) >= threshold;
    HS1 = 'l';
else
    HS1 = 'r';
end

%% flank detection
is_stance = [is_stance(end); is_stance(:)];
is_init_contact = is_stance(2:end) - is_stance(1:end-1);

idx_init_contact = find(is_init_contact == 1);
idx_last_contact = find(is_init_contact == -1);



idx_GC = [idx_init_contact(1):N, 1:idx_init_contact(1)-1]';

















end