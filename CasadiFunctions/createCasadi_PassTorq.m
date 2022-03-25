function [f_PassiveMoments,f_passiveTATorques,f_AllPassiveTorques] = createCasadi_PassTorq(S,model_info)
%% createCasadi_PassTorq.m
%Function to create Casadi functions for passive torques.
%
%INPUT
% - model_info
%
%OUTPUT
%  - Casadi functions
%
%TO DO
% Put hard coded data in Settings (confirm values as well)
% Authors: Ines Vandekerckhove, Dhruv Gupta, KU Leuven
% Date: 30-11-2021 

import casadi.*

%% Passive joint torques (bushing forces/ coodinate limit forces)
K_pass      = SX.sym('K_pass',4);
theta_pass  = SX.sym('theta_pass',2);
qin_pass    = SX.sym('qin_pass',1);
qdotin_pass = SX.sym('qdotin_pass',1);
% theta_pass 1 and 2 are inverted on purpose.
Tau_pass = K_pass(1,1)*exp(K_pass(2,1)*(qin_pass-theta_pass(2,1))) + ...
    K_pass(3,1)*exp(K_pass(4,1)*(qin_pass-theta_pass(1,1))) ...
    - 0.1*qdotin_pass;
f_PassiveMoments = Function('f_PassiveMoments',{K_pass,theta_pass,...
    qin_pass,qdotin_pass},{Tau_pass},{'K_pass','theta_pass',...
    'qin_pass','qdotin_pass'},{'Tau_pass'});

%% Passive torque actuated joint torques (linear stiffnes and damping in joints)
stiff	= SX.sym('stiff',1);
damp	= SX.sym('damp',1);
qin     = SX.sym('qin_pass',1);
qdotin  = SX.sym('qdotin_pass',1);
passTATorques = -stiff * qin - damp * qdotin;
f_passiveTATorques = Function('f_passiveTATorques',{stiff,damp,qin,qdotin}, ...
    {passTATorques},{'stiff','damp','qin','qdotin'},{'passTATorques'});

%% Create function to compute passive moments

Q_SX =  SX.sym('Q_SX',model_info.ExtFunIO.jointi.nq.all,1); % Muscle activations
Qdot_SX = SX.sym('Qdot_SX',model_info.ExtFunIO.jointi.nq.all,1); % Muscle activations

% Get passive joint torques
Tau_passj_muscleActuated = [];
for i=1:model_info.ExtFunIO.jointi.nq.legs_nomtp
    Tau_passj_muscleActuated_var    = f_PassiveMoments(S.k_pass.(model_info.ExtFunIO.jointi.names.legs_nomtp{i}(1:end-2)),...
    S.theta.pass.(model_info.ExtFunIO.jointi.names.legs_nomtp{i}(1:end-2)),...
    Q_SX(model_info.ExtFunIO.jointi.legs_nomtp(i)),...
    Qdot_SX(model_info.ExtFunIO.jointi.legs_nomtp(i)));
    Tau_passj_muscleActuated = [Tau_passj_muscleActuated; Tau_passj_muscleActuated_var];
end
for i=1:model_info.ExtFunIO.jointi.nq.torso
    Tau_passj_muscleActuated_var    = f_PassiveMoments(S.k_pass.(model_info.ExtFunIO.jointi.names.torso{i}),...
    S.theta.pass.(model_info.ExtFunIO.jointi.names.torso{i}),...
    Q_SX(model_info.ExtFunIO.jointi.torso(i)),...
    Qdot_SX(model_info.ExtFunIO.jointi.torso(i)));
    Tau_passj_muscleActuated = [Tau_passj_muscleActuated; Tau_passj_muscleActuated_var];
end

Tau_passj_arms = [];
for i=1:model_info.ExtFunIO.jointi.nq.arms
    Tau_passj_arms_var = f_passiveTATorques(S.stiffnessArm, S.dampingArm, ...
        Q_SX(model_info.ExtFunIO.jointi.armsi(i)), Qdot_SX(model_info.ExtFunIO.jointi.armsi(i)));
    Tau_passj_arms = [Tau_passj_arms; Tau_passj_arms_var];
end

Tau_passj_mtp = [];
for i=1:model_info.ExtFunIO.jointi.nq.mtp
    Tau_passj_mtp_var = f_passiveTATorques(S.kMTP, S.dMTP, ...
    Q_SX(model_info.ExtFunIO.jointi.mtpi(i)), Qdot_SX(model_info.ExtFunIO.jointi.mtpi(i)));
    Tau_passj_mtp = [Tau_passj_mtp; Tau_passj_mtp_var];
end 

Tau_passj_all = [Tau_passj_muscleActuated; Tau_passj_arms; Tau_passj_mtp];

f_AllPassiveTorques = Function('f_AllPassiveTorques',{Q_SX,Qdot_SX}, ...
    {Tau_passj_all},{'Q_SX','Qdot_SX'},{'Tau_passj_all'});

end