function [R] = PostProcess_msk_geometry(model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcess_msk_geometry
%   This function computes the muscle-tendon lengths, velocities and moment
%   arms based on the approximated function. Results are compared against 
%   the OpenSim model.
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

import casadi.*

%% Get approximated geometry from approximated CasADi function
R.muscles.lMT = zeros(N,NMuscle);
R.muscles.vMT = zeros(N,NMuscle);
R.muscles.dM = zeros(N,NMuscle,size(R.kinematics.Qs,2));

for i=1:N

    [lMTj,vMTj,dMj] =  f_casadi.lMT_vMT_dM(R.kinematics.Qs_rad(i,:)',R.kinematics.Qdots_rad(i,:)');

    R.muscles.lMT(i,:) = full(lMTj);
    R.muscles.vMT(i,:) = full(vMTj);
    R.muscles.dM(i,:,:) = full(dMj);

end


%% Get reference geometry from OpenSim model
% [MuscleData] = muscleAnalysisAPI(R.S,model_info.osim_path,model_info,R.kinematics.Qs_rad);
% 
% % Compare
% Delta_lMT = R.muscles.lMT - MuscleData.lMT;
% Delta_dM = R.muscles.dM - MuscleData.dM;
% 
% rmse_lMT = sqrt(mean(Delta_lMT.^2,1));
% rmse_dM = squeeze(sqrt(mean(Delta_dM.^2,1)));
% 
% [il] = find(rmse_lMT > R.S.misc.threshold_lMT_fit);
% [id,jd] = find(squeeze(rmse_dM) > R.S.misc.threshold_dM_fit);
% 
% if ~isempty(il)
%     disp('Issue when reconstructing muscle-tendon length wrt OpenSim model for:')
%     for i=1:length(il)
%         disp(['   ' model_info.muscle_info.muscle_names{il(i)} '(RMSE = ' num2str(rmse_lMT(il(i))*1e3) ' mm)'])
%     end
% end
% 
% if ~isempty(id)
%     disp('Issue when reconstructing muscle-tendon momentarm wrt OpenSim model for:')
%     for i=1:length(id)
%         disp(['   ' model_info.muscle_info.muscle_names{id(i)} ' around ',...
%             model_info.ExtFunIO.coord_names.all{jd(i)} '(RMSE = ' num2str(rmse_dM(id(i),jd(i))*1e3) ' mm)'])
%     end
% end



end