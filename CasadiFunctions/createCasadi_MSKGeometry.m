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
% 
% Last edit by: Menthy Denayer
% Last edit date: 01/September/2026 : Added partial derivative of moment arms to joint angles to f_lMT_vMT_dM function
%
% --------------------------------------------------------------------------
% This file is part of PredSim.
% 
% PredSim: A Framework for Rapid Predictive Simulations of Locomotion
% Copyright (c) 2026 KU Leuven
% 
% PredSim is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Affero General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% 
% PredSim is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public 
% License for more details.
% 
% You should have received a copy of the GNU Affero General Public License 
% along with PredSim. If not, see <https://www.gnu.org/licenses/>.
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
            qin(1,find(muscle_spanning_info_m(i,:)==1)),...
            model_info.muscle_info.polyFit.MuscleInfo.muscle(i).mu);
    end
    
    vMT = jtimes(lMT, qin, qdotin); % v = jacobian(l,q)' * qdot
    dM = - jacobian(lMT, qin);

    % Define casadi function
    if(S.misc.compute_joint_stiffness)
        dMdr = jacobian(dM, qin);                                               % Added by Menthy, to compute the partial derivative of the moment arm to the joint angles, size: (NMuscle*nq) x nq
        f_lMT_vMT_dM = Function('f_lMT_vMT_dM',{qin,qdotin},{lMT,vMT,dM,dMdr}); % Changed by Menthy, added "dMdr" as output for the function, size: (NMuscle*nq) x nq
    else
        f_lMT_vMT_dM = Function('f_lMT_vMT_dM',{qin,qdotin},{lMT,vMT,dM});
    end
    
    % Save function for later use
    f_lMT_vMT_dM.save(msk_geom_path);

end

end
