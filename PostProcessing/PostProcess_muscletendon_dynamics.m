function [R] = PostProcess_muscletendon_dynamics(S,model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcess_muscletendon_dynamics
%   This function computes the muscle-tendon forces, fiber lenghts- and
%   velocities, and tendon lengts- and velocities.
% 
% INPUT:
%   - S -
%   * setting structure S
%
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

N = size(R.Qs,1);
NMuscle = model_info.muscle_info.NMuscle;

tensions = struct_array_to_double_array(model_info.muscle_info.parameters,'specific_tension');

FT = zeros(N,NMuscle);
R.FT = FT;
R.Fce = FT;
R.Fpass = FT;
R.Fiso = FT;
R.lM = FT;
R.lMtilde = FT;
R.vM = FT;
R.vMtilde = FT;
R.lT = FT;
R.vT = FT;

for i=1:N

    [~,FTj,Fcej,Fpassj,Fisoj] = f_casadi.forceEquilibrium_FtildeState_all_tendon(R.a(i,:)',...
        R.FTtilde(i,:)',R.dFTtilde(i,:)',R.lMT(i,:)',R.vMT(i,:)',tensions);

    R.FT(i,:) = full(FTj);
    R.Fce(i,:) = full(Fcej);
    R.Fpass(i,:) = full(Fpassj);
    R.Fiso(i,:) = full(Fisoj);

    [lMj,lMtildej] = f_casadi.FiberLength_TendonForce_tendon(R.FTtilde(i,:)',R.lMT(i,:)');
    
    R.lM(i,:) = full(lMj);
    R.lMtilde(i,:) = full(lMtildej);

    [vMj,vMtildej] = f_casadi.FiberVelocity_TendonForce_tendon(...
        R.FTtilde(i,:)',R.dFTtilde(i,:)',R.lMT(i,:)',R.vMT(i,:)');

    R.vM(i,:) = full(vMj);
    R.vMtilde(i,:) = full(vMtildej);

    [lTj,vTj] = f_casadi.lT_vT(R.FTtilde(i,:)',R.dFTtilde(i,:)',R.lMT(i,:)',R.vMT(i,:)');

    R.lT(i,:) = full(lTj);
    R.vT(i,:) = full(vTj);
end












