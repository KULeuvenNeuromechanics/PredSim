function [f_lMT_vMT_dM] = createCasadi_MSKGeometry(S,model_info)
% --------------------------------------------------------------------------
% createCasadi_MSKGeometry 
%   Function to create Casadi functions for musculoskeletal geometry.
% 
% INPUT:
%   - S -
%   * setting structure S
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% OUTPUT:
%   - f_lMT_vMT_dM -
%   * Casadi function for musculoskeletal geometry.
% 
% Original authors: Lars D'Hondt, Dhruv Gupta, Tom Buurke
% Original date: 01/12/2021
% --------------------------------------------------------------------------


%% Polynomial approximation
import casadi.*

% Check for existing file with polynomial approximation, and load if it
% exists. We only perform muscle analysis and fitting if the result is not 
% yet available, because the analysis takes long.
msk_geom_path = fullfile(S.misc.subject_path,[S.misc.msk_geom_name,'.casadi']);
if isfile(msk_geom_path) && ~S.misc.msk_geom_always_new_fit
    f_lMT_vMT_dM = Function.load(msk_geom_path);
    
elseif strcmpi(S.misc.msk_geom_eq,'polynomials') 
    % Assemble polynomial approximation from coefficients
    muscle_spanning_info_m = model_info.muscle_info.muscle_spanning_joint_info;

    qin     = SX.sym('qin',1,model_info.ExtFunIO.jointi.nq.all);
    qdotin  = SX.sym('qdotin',1,model_info.ExtFunIO.jointi.nq.all);
    lMT     = SX(model_info.muscle_info.NMuscle,1);

    for i=1:size(muscle_spanning_info_m,1)     

        lMT(i,1) = mvpolyval(...
            model_info.muscle_info.polyFit.MuscleInfo.muscle(i).coeff,...
            qin(1,muscle_spanning_info_m(i,:)==1),...
            model_info.muscle_info.polyFit.MuscleInfo.muscle(i).mu);
    end
    
    vMT = jtimes(lMT, qin, qdotin); % v = jacobian(l,q)' * qdot
    dM = - jacobian(lMT, qin);

    % Define casadi function
    f_lMT_vMT_dM = Function('f_lMT_vMT_dM',{qin,qdotin},{lMT,vMT,dM});
    
    % Save function for later use
    f_lMT_vMT_dM.save(msk_geom_path);

end

end
