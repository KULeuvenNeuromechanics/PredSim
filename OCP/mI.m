load('model_info.mat')
% Vectors of indices for later use
residualsi = 1:length(fields(model_info.ExtFunIO.coordi));
ground_pelvisi      = model_info.ExtFunIO.jointi.ground_pelvis; % ground-pelvis
trunki              = model_info.ExtFunIO.jointi.torso; % trunk
legsi               = [model_info.ExtFunIO.jointi.leg_r model_info.ExtFunIO.jointi.leg_l]; % arms
armsi               = [model_info.ExtFunIO.jointi.arm_r model_info.ExtFunIO.jointi.arm_l]; % arms
mtpi                = [model_info.ExtFunIO.jointi.mtp_r model_info.ExtFunIO.jointi.mtp_l]; % mtps
coord_noarmsi       = [ground_pelvisi legsi trunki]; % all but arms
coord_muscleActuated= [legsi trunki]; % all but arms
% Number of degrees of freedom for later use
nq.all      = length(residualsi); % all
nq.abs      = length(ground_pelvisi); % ground-pelvis
nq.trunk    = length(trunki); % trunk
nq.arms     = length(armsi); % arms
nq.mtp      = length(mtpi); % arms
nq.leg      = length(legsi);
nq.muscAct  = length(coord_muscleActuated);

model_info.ExtFunIO.nq = nq;

orderQs = [2*model_info.ExtFunIO.jointi.ground_pelvis(1)-1:2*model_info.ExtFunIO.jointi.ground_pelvis(end),...
    2*model_info.ExtFunIO.jointi.leg_r(1)-1:2*model_info.ExtFunIO.jointi.leg_r(end),...
    2*model_info.ExtFunIO.jointi.leg_l(1)-1:2*model_info.ExtFunIO.jointi.leg_l(end),...
    2*model_info.ExtFunIO.jointi.torso(1)-1:2*model_info.ExtFunIO.jointi.torso(end),...
    2*model_info.ExtFunIO.jointi.arm_r(1)-1:2*model_info.ExtFunIO.jointi.arm_r(end),...
    2*model_info.ExtFunIO.jointi.arm_l(1)-1:2*model_info.ExtFunIO.jointi.arm_l(end)];

orderQsInv = [2*model_info.ExtFunIO.jointi.ground_pelvis(1)-1:2*model_info.ExtFunIO.jointi.ground_pelvis(end),...
    2*model_info.ExtFunIO.jointi.leg_l(1)-1:2*model_info.ExtFunIO.jointi.leg_l(end),...
    2*model_info.ExtFunIO.jointi.leg_r(1)-1:2*model_info.ExtFunIO.jointi.leg_r(end),...
    2*model_info.ExtFunIO.jointi.torso(1)-1:2*model_info.ExtFunIO.jointi.torso(end),...
    2*model_info.ExtFunIO.jointi.arm_l(1)-1:2*model_info.ExtFunIO.jointi.arm_l(end),...
    2*model_info.ExtFunIO.jointi.arm_r(1)-1:2*model_info.ExtFunIO.jointi.arm_r(end)];

orderQsOpp = [2*model_info.ExtFunIO.coordi.pelvis_list-1:2*model_info.ExtFunIO.coordi.pelvis_list,...   
    2*model_info.ExtFunIO.coordi.pelvis_rotation-1:2*model_info.ExtFunIO.coordi.pelvis_rotation,...
    2*model_info.ExtFunIO.coordi.pelvis_tz-1:2*model_info.ExtFunIO.coordi.pelvis_tz,...
    2*model_info.ExtFunIO.coordi.lumbar_bending-1:2*model_info.ExtFunIO.coordi.lumbar_bending,...
    2*model_info.ExtFunIO.coordi.lumbar_rotation-1:2*model_info.ExtFunIO.coordi.lumbar_rotation];


orderArm = [model_info.ExtFunIO.jointi.arm_r,model_info.ExtFunIO.jointi.arm_l]-model_info.ExtFunIO.jointi.arm_r(1)+1;
orderArmInv = [model_info.ExtFunIO.jointi.arm_l,model_info.ExtFunIO.jointi.arm_r]-model_info.ExtFunIO.jointi.arm_r(1)+1;

model_info.ExtFunIO.symQs.orderQs = orderQs;
model_info.ExtFunIO.symQs.orderQsInv = orderQsInv;
model_info.ExtFunIO.symQs.orderQsOpp = orderQsOpp;
model_info.ExtFunIO.symQs.orderArm = orderArm;
model_info.ExtFunIO.symQs.orderArmInv = orderArmInv;

