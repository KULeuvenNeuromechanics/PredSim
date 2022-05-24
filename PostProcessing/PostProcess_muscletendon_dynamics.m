function [R] = PostProcess_muscletendon_dynamics(model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcess_muscletendon_dynamics
%   This function computes the muscle-tendon forces, fiber lenghts- and
%   velocities, and tendon lengts- and velocities.
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
% Original date: 13/May/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

N = size(R.kinematics.Qs,1);
NMuscle = model_info.muscle_info.NMuscle;

tensions = struct_array_to_double_array(model_info.muscle_info.parameters,'specific_tension');

FT = zeros(N,NMuscle);
R.muscles.FT = FT;
R.muscles.Fce = FT;
R.muscles.Fpass = FT;
R.muscles.Fiso = FT;
R.muscles.lM = FT;
R.muscles.lMtilde = FT;
R.muscles.vM = FT;
R.muscles.vMtilde = FT;
R.muscles.lT = FT;
R.muscles.vT = FT;

for i=1:N

    [~,FTj,Fcej,Fpassj,Fisoj] = f_casadi.forceEquilibrium_FtildeState_all_tendon(R.muscles.a(i,:)',...
        R.muscles.FTtilde(i,:)',R.muscles.dFTtilde(i,:)',R.muscles.lMT(i,:)',R.muscles.vMT(i,:)',tensions);

    R.muscles.FT(i,:) = full(FTj);
    R.muscles.Fce(i,:) = full(Fcej);
    R.muscles.Fpass(i,:) = full(Fpassj);
    R.muscles.Fiso(i,:) = full(Fisoj);

    [lMj,lMtildej] = f_casadi.FiberLength_TendonForce_tendon(R.muscles.FTtilde(i,:)',R.muscles.lMT(i,:)');
    
    R.muscles.lM(i,:) = full(lMj);
    R.muscles.lMtilde(i,:) = full(lMtildej);

    [vMj,vMtildej] = f_casadi.FiberVelocity_TendonForce_tendon(...
        R.muscles.FTtilde(i,:)',R.muscles.dFTtilde(i,:)',R.muscles.lMT(i,:)',R.muscles.vMT(i,:)');

    R.muscles.vM(i,:) = full(vMj);
    R.muscles.vMtilde(i,:) = full(vMtildej);

    [lTj,vTj] = f_casadi.lT_vT(R.muscles.FTtilde(i,:)',R.muscles.dFTtilde(i,:)',R.muscles.lMT(i,:)',R.muscles.vMT(i,:)');

    R.muscles.lT(i,:) = full(lTj);
    R.muscles.vT(i,:) = full(vTj);
end












