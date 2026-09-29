function R_left = get_left_gait_cycle(R_right, model_info)
% get_left_gait_cycle
%   Extracts a left gait cycle from PredSim results using GRF-based gait
%   event detection. The function swaps the GRF order to detect the left
%   stance phase, extracts the corresponding kinematics and kinetics, and
%   adjusts the forward position to start at zero.
%
% INPUT:
%   R_right    - PredSim results structure
%   model_info - PredSim model information structure
%
% OUTPUT:
%   R_left     - PredSim results structure containing one left gait cycle

% Original author: Alexandra Duarte Fernandes
% Original date: June 16, 2026

% Last edit by: Ellis Van Can
% Last edit date: September 3, 2026
%% Extract information
dist_trav_opt = R_right.spatiotemp.dist_trav;

GRFs = R_right.ground_reaction;
GRF_right = GRFs.GRF_r;
GRF_left  = GRFs.GRF_l;


%% Detect left gait cycle
% Swap GRF order so the algorithm detects the left stance phase
GRFk_opt = [GRF_left, GRF_right];

[idx_GC, idx_GC_base_forward_offset, ~, ~] = ...
    getStancePhaseSimulation(GRFk_opt, model_info.mass/3);


%% Copy kinematics
R_left.colheaders = R_right.colheaders;

R_left.kinematics.Qs = R_right.kinematics.Qs(idx_GC,:);
R_left.kinematics.Qdots = R_right.kinematics.Qdots(idx_GC,:);
R_left.kinematics.Qddots = R_right.kinematics.Qddots(idx_GC,:);

R_left.kinematics.Qs_rad = R_right.kinematics.Qs_rad(idx_GC,:);
R_left.kinematics.Qdots_rad = R_right.kinematics.Qdots_rad(idx_GC,:);
R_left.kinematics.Qddots_rad = R_right.kinematics.Qddots_rad(idx_GC,:);

%% Copy kinetics
R_left.kinetics.T_ID = R_right.kinetics.T_ID(idx_GC,:);


%% Correct forward position
% Make forward position continuous and start at zero
R_left.kinematics.Qs(idx_GC_base_forward_offset,...
    model_info.ExtFunIO.jointi.base_forward) = ...
    R_left.kinematics.Qs(idx_GC_base_forward_offset,...
    model_info.ExtFunIO.jointi.base_forward) + dist_trav_opt;


R_left.kinematics.Qs(:,model_info.ExtFunIO.jointi.base_forward) = ...
    R_left.kinematics.Qs(:,model_info.ExtFunIO.jointi.base_forward) - ...
    R_left.kinematics.Qs(1,model_info.ExtFunIO.jointi.base_forward);


%% Copy GRFs
R_left.ground_reaction.GRF_r = GRF_right(idx_GC,:);
R_left.ground_reaction.GRF_l = GRF_left(idx_GC,:);


%% Copy spatiotemporal variables
R_left.spatiotemp = R_right.spatiotemp;

end