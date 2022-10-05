function [R] = PostProcess_spatio_temporal(model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcess_spatio_temporal
%   This function computes the spatio-temporal characteristics of the
%   predicted gait.
% 
% INPUT:
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
%   - f_casadi -
%   * Struct containing all casadi functions.
%
%   - R -
%   * struct with simulation results
%
% OUTPUT:
%   - R -
%   * struct with simulation results
% 
% Original author: Lars D'Hondt
% Original date: 19/May/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------



% step length
if isfield(f_casadi,'f_getStepLength')
    [step_length_r,step_length_l] = f_casadi.f_getStepLength(R.kinematics.Qs_rad(1,:),...
        R.kinematics.Qs_rad(end,:));
    R.spatiotemp.step_length_r = full(step_length_r);
    R.spatiotemp.step_length_l = full(step_length_l);
else
    R.spatiotemp.step_length_r = [];
    R.spatiotemp.step_length_l = [];
end

% percentage stance and swing phase
R.spatiotemp.stance_r = length(R.ground_reaction.idx_stance_r)/size(R.kinematics.Qs,1)*100;
R.spatiotemp.swing_r = 100 - R.spatiotemp.stance_r;
R.spatiotemp.stance_l = length(R.ground_reaction.idx_stance_l)/size(R.kinematics.Qs,1)*100;
R.spatiotemp.swing_l = 100 - R.spatiotemp.stance_l;
R.ground_reaction.idx_double_support = intersect(R.ground_reaction.idx_stance_r,...
    R.ground_reaction.idx_stance_l);
R.spatiotemp.double_support = length(R.ground_reaction.idx_double_support)/size(R.kinematics.Qs,1)*100;

% stride frequency
R.spatiotemp.stride_freq = 1/R.time.mesh_GC(end);

% stepwidth
if isfield(R.ground_reaction,'COP_r') && isfield(R.ground_reaction,'COP_l')
    R.spatiotemp.step_width_COP = abs(mean(R.ground_reaction.COP_r(R.ground_reaction.idx_stance_r,3)) ...
        - mean(R.ground_reaction.COP_l(R.ground_reaction.idx_stance_l,3)));
else
    R.spatiotemp.step_width_COP = nan;
end

% distance traveled
R.spatiotemp.dist_trav = R.kinematics.Qs(end,model_info.ExtFunIO.jointi.base_forward) ...
    - R.kinematics.Qs(1,model_info.ExtFunIO.jointi.base_forward);




