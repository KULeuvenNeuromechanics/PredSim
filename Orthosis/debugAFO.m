function [AFO] = debugAFO(init, settings_orthosis)


% create Orthosis object
AFO = Orthosis('AFO',init);

% read settings that were passed from main.m
k_ankle = settings_orthosis.ankle_stiffness; % ankle stiffness in Nm/rad
k_mtp = settings_orthosis.mtp_stiffness; % mtp stiffness in Nm/rad
side = settings_orthosis.left_right; % 'l' for left or 'r' for right

% get joint angles
q_ankle = AFO.var_coord(['ankle_angle_',side]); % ankle angle in rad;
q_mtp = AFO.var_coord(['mtp_angle_',side]); % MTP angle in rad;

if strcmp(side,'l')
    GRF1 = AFO.var_GRF('left_total');
elseif strcmp(side,'r')
    GRF1 = AFO.var_GRF('right_total');
end

GRF2 = AFO.var_GRF(['s1_',side]);

GRF_d = AFO.var_GRF(['s1_',side],'d');

CoM_pos = AFO.var_point('CoM','torso',[-0.1, 0.3, 0]);
CoM_vel = AFO.var_point('CoM','torso',[-0.1, 0.3, 0],'vel');

% calculate moments
T_ankle = [0;0;-k_ankle*q_ankle];
T_mtp = -k_mtp*q_mtp;

% add calculated moments to Orthosis
AFO.addCoordForce(T_mtp, ['mtp_angle_',side])

AFO.addBodyMoment(T_ankle, ['T_exo_shank_',side],['tibia_',side]);
AFO.addBodyMoment(-T_ankle, ['T_exo_foot_',side],['calcn_',side],['tibia_',side]);

end