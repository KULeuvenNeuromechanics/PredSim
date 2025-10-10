function [varargout] = updateSf(usr_height,usr_fingertip_elbow,usr_elbow_shoulder,...
    usr_shoulder_width,usr_hip_knee,usr_knee_ground,usr_foot_length)

%% set defaults for user inputs
default_height = 1.8; % [m]
default_ratio_fingertip_elbow = 1/4;
default_ratio_elbow_shoulder = 1/8;
default_ratio_shoulder_width = 1/4;
default_ratio_hip_knee = 1/4;
default_ratio_knee_ground = 1/4;
default_ratio_foot_length = 1/6;

%% get user input ratios
usr_ratio_fingertip_elbow = usr_fingertip_elbow/usr_height;
usr_ratio_elbow_shoulder = usr_elbow_shoulder/usr_height;
usr_ratio_shoulder_width = usr_shoulder_width/usr_height;
usr_ratio_hip_knee = usr_hip_knee/usr_height;
usr_ratio_knee_ground = usr_knee_ground/usr_height;
usr_ratio_foot_length = usr_foot_length/usr_height;


%% calculate scale factors
sf_foot = usr_ratio_foot_length/default_ratio_foot_length;
sf_low_leg = usr_ratio_knee_ground/default_ratio_knee_ground;
sf_upp_leg = usr_ratio_hip_knee/default_ratio_hip_knee;
sf_upp_arm = usr_ratio_elbow_shoulder/default_ratio_elbow_shoulder;
sf_low_arm = usr_ratio_fingertip_elbow/default_ratio_fingertip_elbow;
sf_shoulder = usr_ratio_shoulder_width/default_ratio_shoulder_width;


%% dimentions for drawing
% foot
dim_foot.length = default_ratio_foot_length*sf_foot;
dim_foot.height = dim_foot.length/4;
dim_foot.width = dim_foot.length/3;
% lower leg
dim_low_leg.length = default_ratio_knee_ground*sf_low_leg - dim_foot.height;
dim_low_leg.width = dim_low_leg.length/3;
% upper leg
dim_upp_leg.length = default_ratio_hip_knee*sf_upp_leg;
dim_upp_leg.width = dim_upp_leg.length/3;
% upper arm
dim_upp_arm.length = default_ratio_elbow_shoulder*sf_upp_arm;
dim_upp_arm.width = dim_upp_arm.length/3;
% hand
dim_hand.length = 1/10;
dim_hand.width = dim_hand.length*0.6;
% lower arm
dim_low_arm.length = (default_ratio_fingertip_elbow*sf_low_arm - dim_hand.length);
dim_low_arm.width = dim_upp_arm.length/3;
% pelvis
dim_pel.width = dim_upp_leg.width*2.2;
dim_pel.height = dim_pel.width*0.4;
% leg length
leg_length = dim_foot.height + dim_low_leg.length + dim_upp_leg.length;
% torso
dim_trs.height = 1-1/6-leg_length;
dim_trs.shwidth = default_ratio_shoulder_width*sf_shoulder;

sf_torso = (1/2) / (1-leg_length);

%% Reurn scale factors if asked
if nargout >= 1
    sf.foot = sf_foot;
    sf.upp_leg = sf_upp_leg;
    sf.low_leg = sf_low_leg;
    sf.torso = sf_torso;
    sf.shoulder = sf_shoulder;
    sf.low_arm = sf_low_arm;
    sf.upp_arm = sf_upp_arm;
    varargout{1} = sf;
end

