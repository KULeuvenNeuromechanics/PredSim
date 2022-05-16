function [R] = PostProcess_msk_geometry(S,model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcess_msk_geometry
%   This function computes the muscle-tendon lengths, velocities and moment
%   arms. Outputs of casadi functions are compared against 
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

lMT = zeros(N,NMuscle);
vMT = lMT;
dM = zeros(N,NMuscle,size(R.Qs,2));

for i=1:N

    [lMTj,vMTj,dMj] =  f_casadi.lMT_vMT_dM(R.Qs_rad(i,:)',R.Qdots_rad(i,:)');

    lMT(i,:) = full(lMTj);
    vMT(i,:) = full(vMTj);
    dM(i,:,:) = full(dMj);

end


% run muscle analysis with predicted motion
[MuscleData] = muscleAnalysisAPI(S,model_info.osim_path,model_info,R.Qs_rad);


%%
Delta_lMT = lMT - MuscleData.lMT;
Delta_dM = dM - MuscleData.dM;

rmse_lMT = sqrt(mean(Delta_lMT.^2,2));
rmse_dM = squeeze(sqrt(mean(Delta_dM.^2,2)));

[il] = find(rmse_lMT >0.003);
[id,jd] = find(squeeze(rmse_dM) >0.003);


%%
R.lMT = lMT;
R.vMT = vMT;
R.dM = dM;

end