function [f_casadi] = createCasadiFunctions(main_path,model_info,N_musc_cross,nq,coord_muscleActuated,armsi,MTparameters_m)
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
% ***Issue: N_musc_cross will be replaced by info from model_info later***
f_casadi = createCasadi_GenHelper(model_info,N_musc_cross);

%% Create Casadi functions for musculoskeletal geometry
f_casadi.lMT_vMT_dM = createCasadi_MSKGeometry(main_path,model_info);

%% Create Casadi functions for muscle contraction dynamics
[f_casadi.forceEquilibrium_FtildeState_all_tendon,f_casadi.FiberLength_TendonForce_tendon...
    ,f_casadi.FiberVelocity_TendonForce_tendon,f_casadi.lT_vT] = createCasadi_ContractDynam(main_path,model_info);

%% Create Casadi functions for passive torques
% *** nq,coord_muscleActuated,armsi to be replaced with info from
% model_info later***
[f_casadi.PassiveMoments,f_casadi.passiveTATorques,f_casadi.AllPassiveTorques]...
    = createCasadi_PassTorq(model_info,nq,coord_muscleActuated,armsi);

%% Create Casadi functions for activation dynamics
[f_casadi.ArmActivationDynamics,f_casadi.TrunkActivationDynamics,f_casadi.MtpActivationDynamics] = createCasadi_ActDynam();

%% Create Casadi functions for metabolic energy.
%*** Replace MTparameters_m by info from model_info later***
[f_casadi.getMetabolicEnergySmooth2004all] = createCasadi_E_Metab(NMuscle,MTparameters_m);
end
