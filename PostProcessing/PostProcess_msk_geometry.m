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

%% Get approximated geometry from CasADi function
R.lMT = zeros(N,NMuscle);
R.vMT = zeros(N,NMuscle);
R.dM = zeros(N,NMuscle,size(R.Qs,2));

for i=1:N

    [lMTj,vMTj,dMj] =  f_casadi.lMT_vMT_dM(R.Qs_rad(i,:)',R.Qdots_rad(i,:)');

    R.lMT(i,:) = full(lMTj);
    R.vMT(i,:) = full(vMTj);
    R.dM(i,:,:) = full(dMj);

end


%% Get reference geometry from OpenSim model
[MuscleData] = muscleAnalysisAPI(S,model_info.osim_path,model_info,R.Qs_rad);

%% Compare
Delta_lMT = R.lMT - MuscleData.lMT;
Delta_dM = R.dM - MuscleData.dM;

rmse_lMT = sqrt(mean(Delta_lMT.^2,1));
rmse_dM = squeeze(sqrt(mean(Delta_dM.^2,1)));

[il] = find(rmse_lMT > S.misc.threshold_lMT_fit);
[id,jd] = find(squeeze(rmse_dM) > S.misc.threshold_dM_fit);

if ~isempty(il)
    disp('Issue when reconstructing muscle-tendon length wrt OpenSim model for:')
    for i=1:length(il)
        disp(['   ' model_info.muscle_info.muscle_names{il(i)} '(RMSE = ' num2str(rmse_lMT(il(i))*1e3) ' mm)'])
    end
end

if ~isempty(id)
    disp('Issue when reconstructing muscle-tendon momentarm wrt OpenSim model for:')
    for i=1:length(id)
        disp(['   ' model_info.muscle_info.muscle_names{id(i)} ' around ',...
            model_info.ExtFunIO.coord_names.all{jd(i)} '(RMSE = ' num2str(rmse_dM(id(i),jd(i))*1e3) ' mm)'])
    end
end

%%

end