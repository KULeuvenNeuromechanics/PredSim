function [f_PassiveStiffnessMoments,f_PassiveDampingMoments,f_LimitTorques,...
    f_AllPassiveTorques,f_AllPassiveTorques_cost] = createCasadi_PassiveMoments(S,model_info)
% --------------------------------------------------------------------------
% createCasadi_PassiveMoments
%   Create CasADi functions to describe the various passive moments
%
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% OUTPUT:
%   - f_PassiveStiffnessMoments -
%   * M = -K*q
%
%   - f_PassiveDampingMoments -
%   * M = -d*qdot
%
%   - f_LimitTorques -
%   * exponential coordinate limit torques
%
%   - f_AllPassiveTorques -
%   * total passive torques
%
%   - f_AllPassiveTorques_cost -
%   * passive torques to use in cost function
%
% Original author: Lars D'Hondt
% Original date: 15/April/2022
%
% Last edit by:
% Last edit date:
% --------------------------------------------------------------------------

import casadi.*

n_coord = model_info.ExtFunIO.jointi.nq.all;
passive_moment_info = model_info.passive_moment_info;

%% Coordinate limit torque 
K_pass      = SX.sym('K_pass',4);
theta_pass  = SX.sym('theta_pass',2);
qin_pass    = SX.sym('qin_pass',1);

% theta_pass 1 and 2 are inverted on purpose.
Tau_pass = K_pass(1,1)*exp(K_pass(2,1)*(qin_pass-theta_pass(2,1))) + ...
    K_pass(3,1)*exp(K_pass(4,1)*(qin_pass-theta_pass(1,1)));

f_limit_torque = Function('f_limit_torque',{K_pass,theta_pass,...
    qin_pass},{Tau_pass},{'K_pass','theta_pass',...
    'qin_pass'},{'Tau_pass'});

%% all coordinates

q = SX.sym('q',n_coord);
qdot = SX.sym('qdot',n_coord);
tau_k = SX(n_coord,1);
tau_d = SX(n_coord,1);
tau_lim = SX(n_coord,1);

for i=1:n_coord
    if passive_moment_info(i).stiffness_coeff~=0
        tau_k(i) = -passive_moment_info(i).stiffness_coeff*q(i);
    end
    if passive_moment_info(i).damping_coeff~=0
        tau_d(i) = -passive_moment_info(i).damping_coeff*qdot(i);
    end

    if ~isempty(passive_moment_info(i).K_pass) && ~isempty(passive_moment_info(i).theta_pass)
        tau_lim(i) = f_limit_torque(passive_moment_info(i).K_pass,passive_moment_info(i).theta_pass,q(i));

    end
end

if S.weights.pass_torq_includes_damping
    tau_J = tau_lim + tau_d;
else
    tau_J = tau_lim;
end
tau_J1 = tau_J(model_info.ExtFunIO.jointi.limitTorque);
tau_tot = tau_k + tau_d + tau_lim;

%% create functions
f_PassiveStiffnessMoments = Function('f_PassiveStiffnessMoments',{q},{tau_k},{'q'},{'tau_k'});

f_PassiveDampingMoments = Function('f_PassiveDampingMoments',{qdot},{tau_d},{'qdot'},{'tau_d'});

f_LimitTorques = Function('f_LimitTorques',{q},{tau_lim},{'q'},{'tau_lim'});

f_AllPassiveTorques = Function('f_AllPassiveTorques',{q,qdot},{tau_tot},{'q','qdot'},{'tau_tot'});

f_AllPassiveTorques_cost = Function('f_AllPassiveTorques_cost',{q,qdot},{tau_J1},{'q','qdot'},{'tau_J'});





















