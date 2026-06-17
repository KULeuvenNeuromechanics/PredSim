function [R] = PostProcess_ligament_geometry(model_info, f_casadi, R)

% --------------------------------------------------------------------------
% PostProcess_ligament
%   This function computes the ligament parameters from - based on
%   PostProcessing_msk_geometry() ;
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
%
% --------------------------------------------------------------------------


N = size(R.kinematics.Qs,1);
NLigs = model_info.ligament_info.NLigament;

import casadi.*

%% Get approximated geometry from approximated CasADi function
R.ligaments.lLig = zeros(N,NLigs);
R.ligaments.vLig = zeros(N,NLigs);
R.ligaments.mLig = zeros(N,size(R.kinematics.Qs,2));
R.ligaments.fLig = zeros(N,NLigs);


for i=1:N

    [lLig,vLig,fLig] =  f_casadi.ligamentLengthForce(R.kinematics.Qs_rad(i,:)',R.kinematics.Qdots_rad(i,:)');

    R.ligaments.lLig(i,:) = full(lLig);
    R.ligaments.vLig(i,:) = full(vLig);
    R.ligaments.fLig(i,:) = full(fLig);
    
    [mLig] = f_casadi.ligamentMoment(R.kinematics.Qs_rad(i,:)');
    R.ligaments.mLig(i,:) = full(mLig);

end













end

