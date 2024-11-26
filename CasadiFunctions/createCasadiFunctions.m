function [f_casadi] = createCasadiFunctions(S,model_info)
% --------------------------------------------------------------------------
%createCasadiFunctions.m
%   Overview function from which al casadi functions are created
% 
% INPUT:
%   - S -
%   * setting structure S
%  
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% OUTPUT:
%   - f_casadi -
%   * Struct containing all casadi functions.
% 
% Original author: Tom Buurke
% Original date: 02/12/2021
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

%% Create generic casadi functions
f_casadi = createCasadi_GenHelper(S,model_info);

%% Create Casadi functions for musculoskeletal geometry
f_casadi.lMT_vMT_dM = createCasadi_MSKGeometry(S,model_info);

%% Create Casadi functions for muscle contraction dynamics
[forceEquilibrium_FtildeState_all_tendon, FiberLength_TendonForce_tendon,...
    FiberVelocity_TendonForce_tendon,lT_vT] = createCasadi_ContractDynam(S,model_info);

f_casadi.forceEquilibrium_FtildeState_all_tendon = forceEquilibrium_FtildeState_all_tendon;
f_casadi.FiberLength_TendonForce_tendon = FiberLength_TendonForce_tendon;
f_casadi.FiberVelocity_TendonForce_tendon = FiberVelocity_TendonForce_tendon;
f_casadi.lT_vT = lT_vT;

%% Create Casadi functions for angle-moment relation from ligaments
[f_casadi.ligamentMoment, f_casadi.ligamentMoment_single,...
    f_casadi.ligamentMoment_multi, f_casadi.ligamentLengthForce] =...
    createCasadi_Ligaments(S,model_info);

%% Create Casadi functions for passive torques
[f_casadi.PassiveStiffnessMoments, f_casadi.PassiveDampingMoments,...
    f_casadi.LimitTorques, f_casadi.AllPassiveTorques,...
    f_casadi.AllPassiveTorques_cost] = createCasadi_PassiveMoments(S,model_info);

%% Create Casadi functions for activation dynamics
if model_info.ExtFunIO.jointi.nq.torqAct > 0
    [f_casadi.ActuatorActivationDynamics] = createCasadi_ActDynam(S,model_info);
end

%% Create Casadi functions for metabolic energy.
[f_casadi.getMetabolicEnergySmooth2004all] = createCasadi_E_Metab(S,model_info);

%% Create Casadi function to get step length
if ~isempty(S.bounds.SLL.upper) || ~isempty(S.bounds.SLR.upper)
    [f_casadi.f_getCalcnOriginInWorldFrame,f_casadi.f_getStepLength] = ...
        createCasadi_StepLength(S,model_info);
end

%% Create Casadi function for Synergy weight tracking term
if (S.subject.synergies) && (S.subject.TrackSynW)
    [f_casadi.TrackSynW] = createCasadi_TrackSynW(S,model_info);
end
%% Create Casadi functions for orthoses
[f_casadi.f_orthosis_mesh_k, f_casadi.f_orthosis_mesh_all,...
    f_casadi.separate_orthoses] = createCasadi_orthosis(S,model_info);

end