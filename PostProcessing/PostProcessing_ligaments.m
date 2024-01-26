function [R] = PostProcessing_ligaments(model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcessing_ligaments
%   Calculate lengths, forces, powers of ligaments
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
% Original date: (17/April/2023
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

if model_info.ligament_info.NLigament == 0
    % there are no ligaments
    return
end

R.colheaders.ligaments = model_info.ligament_info.ligament_names;

[L, v, F]= f_casadi.ligamentLengthForce(R.kinematics.Qs_rad', R.kinematics.Qdots_rad');

R.ligaments.l = full(L');
R.ligaments.v = full(v');
R.ligaments.F = full(F');

CSA = struct_array_to_double_array(model_info.ligament_info.parameters,'cross_section_area');
ls = struct_array_to_double_array(model_info.ligament_info.parameters,'slack_length');

R.ligaments.strain = (R.ligaments.l./ls'-1)*100;
R.ligaments.stress = R.ligaments.F./CSA';

R.ligaments.moment = full(f_casadi.ligamentMoment(R.kinematics.Qs_rad'))';
R.ligaments.power = R.ligaments.moment.*R.kinematics.Qdots_rad;


end % end of function