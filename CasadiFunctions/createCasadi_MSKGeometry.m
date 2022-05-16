function [f_lMT_vMT_dM] = createCasadi_MSKGeometry(S,model_info)
%createCasadi_MSKGeometry 
%   Function to create Casadi functions for musculoskeletal geometry.
% 
% INPUT:
%   MainPath
%   * Main path
%  
%   model_info
%   * Model info struct
% 
% OUTPUT:
%   f_lMT_vMT_dM
%   * Casadi function for musculoskeletal geometry.
% 
% Original authors: Lars D'Hondt, Dhruv Gupta, Tom Buurke
% Original date: 01/12/2021


%% Polynomial approximation
import casadi.*

% Check for existing file with polynomial approximation, and load if it
% exists. We only perform muscle analysis and fitting if the result is not 
% yet available, because the analysis takes long.
if isfile(fullfile(S.misc.subject_path,S.misc.msk_geom_name))
    f_lMT_vMT_dM = Function.load(fullfile(S.misc.subject_path,S.misc.msk_geom_name));
else
    % Assemble polynomial approximation from coefficients
    muscle_spanning_info_m = model_info.muscle_info.muscle_spanning_joint_info(:,:);
    MuscleInfo_m.muscle    = model_info.muscle_info.polyFit.MuscleInfo.muscle(:);
    NMuscle_pol = size(muscle_spanning_info_m,1);
    NJoint_pol = size(muscle_spanning_info_m,2);
    qin     = SX.sym('qin',1,NJoint_pol);
    qdotin  = SX.sym('qdotin',1,NJoint_pol);
    lMT     = SX(NMuscle_pol,1);
    vMT     = SX(NMuscle_pol,1);
    dM      = SX(NMuscle_pol,NJoint_pol);
    load nCoeffMat
    load expoVal_all
    for i=1:NMuscle_pol      
        index_dof_crossing  = find(muscle_spanning_info_m(i,:)==1);
        nr_dof_crossing     = length(index_dof_crossing); 
        order               = MuscleInfo_m.muscle(i).order;
        [mat,diff_mat_q]    = n_art_mat_3_cas_SX_7(qin(1,index_dof_crossing),...
            nCoeffMat(order,nr_dof_crossing),expoVal_all{order,nr_dof_crossing});
        lMT(i,1)            = mat'*MuscleInfo_m.muscle(i).coeff;
        vMT(i,1)            = 0;
        dM(i,1:NJoint_pol)  = 0;
        for dof_nr = 1:nr_dof_crossing
            dM(i,index_dof_crossing(dof_nr)) = ...
                (-(diff_mat_q(:,dof_nr)))'*MuscleInfo_m.muscle(i).coeff;
            vMT(i,1) = vMT(i,1) + (-dM(i,index_dof_crossing(dof_nr))*...
                qdotin(1,index_dof_crossing(dof_nr)));
        end 
    end
    % Define casadi function
    f_lMT_vMT_dM = Function('f_lMT_vMT_dM',{qin,qdotin},{lMT,vMT,dM});
    
    % Save function for later use
    f_lMT_vMT_dM.save(fullfile(S.misc.subject_path,S.misc.msk_geom_name));
end

end