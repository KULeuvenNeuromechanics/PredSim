function [R] = PostProcess_passive_moments(S,model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcess_passive_moments
%   This function computes the different contributions to the passive joint
%   moments.
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
% Original date: 19/May/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


N = size(R.kinematics.Qs,1);

R.kinetics.T_spring = zeros(N,model_info.ExtFunIO.jointi.nq.all);
R.kinetics.T_damping = zeros(N,model_info.ExtFunIO.jointi.nq.all);
R.kinetics.T_limit = zeros(N,model_info.ExtFunIO.jointi.nq.all);

for i=1:N
    Tki = f_casadi.PassiveStiffnessMoments(R.kinematics.Qs_rad(i,:)');
    R.kinetics.T_spring(i,:) = full(Tki)';

    Tdi = f_casadi.PassiveDampingMoments(R.kinematics.Qdots_rad(i,:)');
    R.kinetics.T_damping(i,:) = full(Tdi)';

    Tli = f_casadi.LimitTorques(R.kinematics.Qs_rad(i,:)');
    R.kinetics.T_limit(i,:) = full(Tli)';

end




