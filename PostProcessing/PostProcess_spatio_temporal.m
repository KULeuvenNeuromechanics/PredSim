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

%%% Tom edit STAsym
%Find right heelstrike
dummy_sl_idx_R = find(diff(find(R.ground_reaction.GRF_r(:,2)>threshold))>1);
if isempty(dummy_sl_idx_R)
    temp_sl_idx_R = find(R.ground_reaction.GRF_r(:,2)>threshold,1,'first');
elseif length(dummy_sl_idx_R)>1 % Fix for multiple foot ground-contacts
    disp('Warning: Multiple ground-contacts!')
    [~,HSRMax]=max(R.ground_reaction.GRF_r(dummy_sl_idx_R,2));
    temp_sl_idx_R = find(R.ground_reaction.GRF_r(:,2)>threshold,dummy_sl_idx_R(HSRMax)+1);
else
    temp_sl_idx_R = find(R.ground_reaction.GRF_r(:,2)>threshold,dummy_sl_idx_R+1);
end
HSR = temp_sl_idx_R(end);

%Find left heelstrike
threshold=20;
dummy_sl_idx_L = find(diff(find(R.ground_reaction.GRF_l(:,2)>threshold))>1);
if isempty(dummy_sl_idx_L)
    temp_sl_idx_L = find(R.ground_reaction.GRF_l(:,2)>threshold,1,'first');
elseif length(dummy_sl_idx_L)>1 % Fix for multiple foot ground-contacts
    disp('Warning: Multiple ground-contacts!')
    [~,HSLMax]=max(R.ground_reaction.GRF_l(dummy_sl_idx_L,2));
    temp_sl_idx_L = find(R.ground_reaction.GRF_l(:,2)>threshold,dummy_sl_idx_L(HSLMax)+1);
else
    temp_sl_idx_L = find(R.ground_reaction.GRF_l(:,2)>threshold,dummy_sl_idx_L+1);
end
HSL = temp_sl_idx_L(end);

%Step time asymmetry based on heelstrike
if HSL > HSR %Right step occurs first
    R.spatiotemp.steptime_r = HSL-HSR;
    R.spatiotemp.steptime_l = R.S.solver.N_meshes-R.spatiotemp.steptime_r;
elseif HSR > HSL %Left step occurs first
    HSL=1;
    R.spatiotemp.steptime_l = HSR-HSL;
    R.spatiotemp.steptime_r = R.S.solver.N_meshes-R.spatiotemp.steptime_l;
end
R.spatiotemp.steptime_asym = (R.spatiotemp.steptime_l - R.spatiotemp.steptime_r) / (R.spatiotemp.steptime_l + R.spatiotemp.steptime_r);
R.spatiotemp.steptimesec_l = R.spatiotemp.steptime_l/100*R.time.mesh(end);
R.spatiotemp.steptimesec_r = R.spatiotemp.steptime_r/100*R.time.mesh(end);
R.spatiotemp.steptimesec_asym = (R.spatiotemp.steptimesec_l - R.spatiotemp.steptimesec_r) / (R.spatiotemp.steptimesec_l + R.spatiotemp.steptimesec_r);
%%% Tom edit STAsym end

% stride frequency
R.spatiotemp.stride_freq = 1/R.time.mesh_GC(end);

% stepwidth
R.spatiotemp.step_width_COP = abs(mean(R.ground_reaction.COP_r(R.ground_reaction.idx_stance_r,3)) ...
    - mean(R.ground_reaction.COP_l(R.ground_reaction.idx_stance_l,3)));

% distance traveled
R.spatiotemp.dist_trav = R.kinematics.Qs(end,model_info.ExtFunIO.jointi.base_forward) ...
    - R.kinematics.Qs(1,model_info.ExtFunIO.jointi.base_forward);




