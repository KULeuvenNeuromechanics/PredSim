function S = updateS_passiveJointTorqueProperties(S)
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

S.stiffnessArm = 0; % CONFIRM VALUE
S.dampingArm = 0.1; % CONFIRM VALUE

S.kMTP = 1.5/(pi/180)/5; % CONFIRM VALUE
S.dMTP = 0.5; % CONFIRM VALUE
