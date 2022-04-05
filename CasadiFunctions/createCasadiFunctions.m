function [f_casadi] = createCasadiFunctions(S,model_info)
%createCasadiFunctions.m
%   Overview function from which al casadi functions are created
% 
% INPUT:
%   main_path
%   * Main path
%  
%   model_info
%   * Struct containing all model_info
% 
% OUTPUT:
%   f_casadi
%   * Struct containing all casadi functions.
% 
% Original author: Tom Buurke
% Original date: 02/12/2021

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

%% Create Casadi functions for passive torques
[f_casadi.PassiveMoments,f_casadi.passiveTATorques,...
    f_casadi.AllPassiveTorques] = createCasadi_PassTorq(S,model_info);

%% Create Casadi functions for activation dynamics
[f_casadi.ArmActivationDynamics,f_casadi.TrunkActivationDynamics,...
    f_casadi.MtpActivationDynamics] = createCasadi_ActDynam(S,model_info);

%% Create Casadi functions for metabolic energy.
[f_casadi.getMetabolicEnergySmooth2004all] = createCasadi_E_Metab(S,model_info);
end