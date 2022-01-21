function [f_PassiveMoments,f_passiveTATorques,f_AllPassiveTorques] = createCasadi_PassTorq(model_info,S)
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

%% Passive joint torques
% We extract the parameters for the passive torques of the lower limbs and
% the trunk

% NEED TO CHANGE THEM IN THE BEGINNING, AND CALL the muscle actuated joints as ExtFunIO.jointi.names.muscleActuated
S.k_pass.hip_flexion = [-2.44 5.05 1.51 -21.88]';
S.theta.pass.hip_flexion = [-0.6981 1.81]';
S.k_pass.hip_adduction = [-0.03 14.94 0.03 -14.94]';
S.theta.pass.hip_adduction = [-0.5 0.5]';
S.k_pass.hip_rotation = [-0.03 14.94 0.03 -14.94]';
S.theta.pass.hip_rotation = [-0.92 0.92]';
S.k_pass.knee_angle = [-6.09 33.94 11.03 -11.33]';
S.theta.pass.knee_angle = [-2.4 0.13]';
S.k_pass.ankle_angle = [-2.03 38.11 0.18 -12.12]';
S.theta.pass.ankle_angle = [-0.74 0.52]';
S.k_pass.subtalar_angle = [-60.21 16.32 60.21 -16.32]';
S.theta.pass.subtalar_angle = [-0.65 0.65]';
S.k_pass.mtp_angle = [-0.9 14.87 0.18 -70.08]';
S.theta.pass.mtp_angle = [0 65/180*pi]';
S.k_pass.lumbar_extension = [-0.35 30.72 0.25 -20.36]';
S.theta.pass.lumbar_extension = [-0.5235987755982988 0.17]';
S.k_pass.lumbar_bending = [-0.25 20.36 0.25 -20.36]';
S.theta.pass.lumbar_bending = [-0.3490658503988659 0.3490658503988659]';
S.k_pass.lumbar_rotation = [-0.25 20.36 0.25 -20.36]';
S.theta.pass.lumbar_rotation = [-0.3490658503988659 0.3490658503988659]';

%% Create function to compute passive moments

Q_SX =  SX.sym('Q_SX',model_info.ExtFunIO.nq.all,1); % Muscle activations
Qdot_SX = SX.sym('Q_SX',model_info.ExtFunIO.nq.all,1); % Muscle activations

%%% Add for loop or adapt computations as a nmatrix rather than individual
%%% coordinate

% Get passive joint torques

Tau_passj_muscleActuated = [];
for i=1:model_info.ExtFunIO.nq.muscleActuated

%   coord_muscleActuated(i)
    Tau_passj_muscleActuated_var    = f_PassiveMoments(S.k_pass.(model_info.ExtFunIO.jointi.names.muscleActuated{i}),...
    S.theta.pass.(model_info.ExtFunIO.jointi.names.muscleActuated{i}),Q_SX(model_info.ExtFunIO.jointi.muscleActuated(i)),...
    Qdot_SX(model_info.ExtFunIO.jointi.muscleActuated(i)));
    Tau_passj_muscleActuated = [Tau_passj_muscleActuated; Tau_passj_muscleActuated_var];

% Tau_passj.hip.flex.l    = f_PassiveMoments(k_pass.hip.flex,...
%     theta.pass.hip.flex,Q_SX(jointi.hip_flex.l),...
%     Qdot_SX(jointi.hip_flex.l));
% Tau_passj.hip.flex.r    = f_PassiveMoments(k_pass.hip.flex,...
%     theta.pass.hip.flex,Q_SX(jointi.hip_flex.r),...
%     Qdot_SX(jointi.hip_flex.r));
% Tau_passj.hip.add.l     = f_PassiveMoments(k_pass.hip.add,...
%     theta.pass.hip.add,Q_SX(jointi.hip_add.l),...
%     Qdot_SX(jointi.hip_add.l));
% Tau_passj.hip.add.r     = f_PassiveMoments(k_pass.hip.add,...
%     theta.pass.hip.add,Q_SX(jointi.hip_add.r),...
%     Qdot_SX(jointi.hip_add.r));
% Tau_passj.hip.rot.l     = f_PassiveMoments(k_pass.hip.rot,...
%     theta.pass.hip.rot,Q_SX(jointi.hip_rot.l),...
%     Qdot_SX(jointi.hip_rot.l));
% Tau_passj.hip.rot.r     = f_PassiveMoments(k_pass.hip.rot,...
%     theta.pass.hip.rot,Q_SX(jointi.hip_rot.r),...
%     Qdot_SX(jointi.hip_rot.r));
% Tau_passj.knee.l        = f_PassiveMoments(k_pass.knee,...
%     theta.pass.knee,Q_SX(jointi.knee.l),...
%     Qdot_SX(jointi.knee.l));
% Tau_passj.knee.r        = f_PassiveMoments(k_pass.knee,...
%     theta.pass.knee,Q_SX(jointi.knee.r),...
%     Qdot_SX(jointi.knee.r));
% Tau_passj.ankle.l       = f_PassiveMoments(k_pass.ankle,...
%     theta.pass.ankle,Q_SX(jointi.ankle.l),...
%     Qdot_SX(jointi.ankle.l));
% Tau_passj.ankle.r       = f_PassiveMoments(k_pass.ankle,...
%     theta.pass.ankle,Q_SX(jointi.ankle.r),...
%     Qdot_SX(jointi.ankle.r));
% Tau_passj.subt.l       = f_PassiveMoments(k_pass.subt,...
%     theta.pass.subt,Q_SX(jointi.subt.l),...
%     Qdot_SX(jointi.subt.l));
% Tau_passj.subt.r       = f_PassiveMoments(k_pass.subt,...
%     theta.pass.subt,Q_SX(jointi.subt.r),...
%     Qdot_SX(jointi.subt.r));
% Tau_passj.trunk.ext     = f_PassiveMoments(k_pass.trunk.ext,...
%     theta.pass.trunk.ext,Q_SX(jointi.trunk.ext),...
%     Qdot_SX(jointi.trunk.ext));
% Tau_passj.trunk.ben     = f_PassiveMoments(k_pass.trunk.ben,...
%     theta.pass.trunk.ben,Q_SX(jointi.trunk.ben),...
%     Qdot_SX(jointi.trunk.ben));
% Tau_passj.trunk.rot     = f_PassiveMoments(k_pass.trunk.rot,...
%     theta.pass.trunk.rot,Q_SX(jointi.trunk.rot),...
%     Qdot_SX(jointi.trunk.rot));
end

% NEED TO ADD ARM STIFFNESS AND DAMPING IN THE BEGINNING
S.stiffnessArm = 0; % CONFIRM VALUE
S.dampingArm = 0.1; % CONFIRM VALUE
Tau_passj_arms = [];

for i=1:model_info.ExtFunIO.nq.arms
    
Tau_passj_arms_var = f_passiveTATorques(S.stiffnessArm, S.dampingArm, ...
    Q_SX(model_info.ExtFunIO.jointi.armsi(i)), Qdot_SX(model_info.ExtFunIO.jointi.armsi(i)));
Tau_passj_arms = [Tau_passj_arms; Tau_passj_arms_var];

% Tau_passj.sh_flex.l = f_passiveTATorques(stiffnessArm, dampingArm, ...
%     Q_SX(jointi.sh_flex.l), Qdot_SX(jointi.sh_flex.l));
% Tau_passj.sh_add.l = f_passiveTATorques(stiffnessArm, dampingArm, ...
%     Q_SX(jointi.sh_add.l), Qdot_SX(jointi.sh_add.l));
% Tau_passj.sh_rot.l = f_passiveTATorques(stiffnessArm, dampingArm, ...
%     Q_SX(jointi.sh_rot.l), Qdot_SX(jointi.sh_rot.l));
% Tau_passj.sh_flex.r = f_passiveTATorques(stiffnessArm, dampingArm, ...
%     Q_SX(jointi.sh_flex.r), Qdot_SX(jointi.sh_flex.r));
% Tau_passj.sh_add.r = f_passiveTATorques(stiffnessArm, dampingArm, ...
%     Q_SX(jointi.sh_add.r), Qdot_SX(jointi.sh_add.r));
% Tau_passj.sh_rot.r = f_passiveTATorques(stiffnessArm, dampingArm, ...
%     Q_SX(jointi.sh_rot.r), Qdot_SX(jointi.sh_rot.r));
% Tau_passj.elb.l = f_passiveTATorques(stiffnessArm, dampingArm, ...
%     Q_SX(jointi.elb.l), Qdot_SX(jointi.elb.l));
% Tau_passj.elb.r = f_passiveTATorques(stiffnessArm, dampingArm, ...
%     Q_SX(jointi.elb.r), Qdot_SX(jointi.elb.r));
% Tau_passj.arm = [Tau_passj.sh_flex.l, Tau_passj.sh_add.l, ...
%     Tau_passj.sh_rot.l, Tau_passj.sh_flex.r, Tau_passj.sh_add.r, ...
%     Tau_passj.sh_rot.r, Tau_passj.elb.l, Tau_passj.elb.r];
end

% NEED TO ADD MTP STIFFNESS AND DAMPING IN THE BEGINNING
S.kMTP = 1.5/(pi/180)/5; % CONFIRM VALUE
S.dMTP = 0.5; % CONFIRM VALUE

Tau_passj_mtp = [];
for i=model_info.ExtFunIO.nq.mtp
%    mtpi(i)
    Tau_passj_mtp_var = f_passiveTATorques(S.kMTP, S.dMTP, ...
    Q_SX(model_info.ExtFunIO.jointi.mtpi(i)), Qdot_SX(model_info.ExtFunIO.jointi.mtpi(i)));
    Tau_passj_mtp = [Tau_passj_mtp; Tau_passj_mtp_var];

% Tau_passj.mtp.l = f_passiveTATorques(Settings.kMTP, Settings.dMTP, ...
%     Q_SX(jointi.mtp.l), Qdot_SX(jointi.mtp.l));
% Tau_passj.mtp.r = f_passiveTATorques(Settings.kMTP, Settings.dMTP, ...
%     Q_SX(jointi.mtp.r), Qdot_SX(jointi.mtp.r));
% Tau_passj.mtp.all = [Tau_passj.mtp.l, Tau_passj.mtp.r];
end 

Tau_passj_all = [Tau_passj_muscleActuated, Tau_passj_arms, Tau_passj_mtp];
% Tau_passj_all = [Tau_passj.hip.flex.l,Tau_passj.hip.flex.r,...
%     Tau_passj.hip.add.l,Tau_passj.hip.add.r,...
%     Tau_passj.hip.rot.l,Tau_passj.hip.rot.r,...
%     Tau_passj.knee.l,Tau_passj.knee.r,Tau_passj.ankle.l,...
%     Tau_passj.ankle.r,Tau_passj.subt.l,Tau_passj.subt.r,...
%     Tau_passj.mtp.all,Tau_passj.trunk.ext,Tau_passj.trunk.ben,...
%     Tau_passj.trunk.rot,Tau_passj.arm]';

f_AllPassiveTorques = Function('f_AllPassiveTorques',{Q_SX,Qdot_SX}, ...
    {Tau_passj_all},{'Q_SX','Qdot_SX'},{'Tau_passj_all'});

end