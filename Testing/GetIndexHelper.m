function [model_info] = GetIndexHelper(S,model_info)
% TO DO
% Why different k for calf muslces?

residualsi = 1:length(fields(model_info.ExtFunIO.coordi));
model_info.ExtFunIO.jointi.legsi = [model_info.ExtFunIO.jointi.leg_r model_info.ExtFunIO.jointi.leg_l];
model_info.ExtFunIO.jointi.armsi = [model_info.ExtFunIO.jointi.arm_r model_info.ExtFunIO.jointi.arm_l];
model_info.ExtFunIO.jointi.mtpi  = [model_info.ExtFunIO.jointi.mtp_r model_info.ExtFunIO.jointi.mtp_l];
model_info.ExtFunIO.jointi.legs_nomtp = setdiff(model_info.ExtFunIO.jointi.legsi,model_info.ExtFunIO.jointi.mtpi);
model_info.ExtFunIO.jointi.noarmsi = [model_info.ExtFunIO.jointi.ground_pelvis model_info.ExtFunIO.jointi.legsi model_info.ExtFunIO.jointi.torso]; % all but arms
model_info.ExtFunIO.jointi.legs_torso = [model_info.ExtFunIO.jointi.legsi model_info.ExtFunIO.jointi.torso];
model_info.ExtFunIO.jointi.muscleActuated = setdiff(model_info.ExtFunIO.jointi.legs_torso,model_info.ExtFunIO.jointi.mtpi);

% Number of degrees of freedom for later use
nq.all      = length(residualsi); % all
nq.abs      = length(model_info.ExtFunIO.jointi.ground_pelvis);
nq.torso    = length(model_info.ExtFunIO.jointi.torso); % trunk
nq.arms     = length(model_info.ExtFunIO.jointi.armsi); % arms
nq.mtp      = length(model_info.ExtFunIO.jointi.mtpi);
nq.leg      = length(model_info.ExtFunIO.jointi.legsi);
nq.legs_torso = length(model_info.ExtFunIO.jointi.legs_torso);
nq.noarms   = length(model_info.ExtFunIO.jointi.noarmsi);
nq.muscleActuated = length(model_info.ExtFunIO.jointi.muscleActuated);
nq.legs_nomtp = length(model_info.ExtFunIO.jointi.legs_nomtp);

model_info.ExtFunIO.nq = nq;

model_info = addCoordNames(model_info,'leg_r');
model_info = addCoordNames(model_info,'leg_l');
model_info = addCoordNames(model_info,'arm_r');
model_info = addCoordNames(model_info,'arm_l');
model_info = addCoordNames(model_info,'ground_pelvis');
model_info = addCoordNames(model_info,'torso');
model_info = addCoordNames(model_info,'muscleActuated');
model_info = addCoordNames(model_info,'legs_torso');
model_info = addCoordNames(model_info,'legsi');
model_info = addCoordNames(model_info,'armsi');
model_info = addCoordNames(model_info,'noarmsi');
model_info = addCoordNames(model_info,'mtpi');
model_info = addCoordNames(model_info,'legs_nomtp');

orderQs = 1:2*nq.all;

coords = fields(model_info.ExtFunIO.coordi);
for i=1:nq.all
    if strcmp(coords{i}(end-1:end),'_r')
        coordInv{i} = [coords{i}(1:end-2) '_l'];
        orderInv(i) = model_info.ExtFunIO.coordi.(coordInv{i});
    elseif strcmp(coords{i}(end-1:end),'_l')
        coordInv{i} = [coords{i}(1:end-2) '_r'];
        orderInv(i) = model_info.ExtFunIO.coordi.(coordInv{i});
    else
        orderInv(i) = model_info.ExtFunIO.coordi.(coords{i});
    end
end
orderQsInv(1:2:2*nq.all) = (2*orderInv)-1;
orderQsInv(2:2:2*nq.all) = 2*orderInv;
    
% orderQsInv = [2*model_info.ExtFunIO.jointi.ground_pelvis(1)-1:2*model_info.ExtFunIO.jointi.ground_pelvis(end),...
%     2*model_info.ExtFunIO.jointi.leg_l(1)-1:2*model_info.ExtFunIO.jointi.leg_l(end),...
%     2*model_info.ExtFunIO.jointi.leg_r(1)-1:2*model_info.ExtFunIO.jointi.leg_r(end),...
%     2*model_info.ExtFunIO.jointi.torso(1)-1:2*model_info.ExtFunIO.jointi.torso(end),...
%     2*model_info.ExtFunIO.jointi.arm_l(1)-1:2*model_info.ExtFunIO.jointi.arm_l(end),...
%     2*model_info.ExtFunIO.jointi.arm_r(1)-1:2*model_info.ExtFunIO.jointi.arm_r(end)];

orderQsOpp1 = [2*model_info.ExtFunIO.coordi.pelvis_list-1:2*model_info.ExtFunIO.coordi.pelvis_list,...   
    2*model_info.ExtFunIO.coordi.pelvis_rotation-1:2*model_info.ExtFunIO.coordi.pelvis_rotation,...
    2*model_info.ExtFunIO.coordi.pelvis_tz-1:2*model_info.ExtFunIO.coordi.pelvis_tz,...
    2*model_info.ExtFunIO.coordi.lumbar_bending-1:2*model_info.ExtFunIO.coordi.lumbar_bending,...
    2*model_info.ExtFunIO.coordi.lumbar_rotation-1:2*model_info.ExtFunIO.coordi.lumbar_rotation];

orderArm = [model_info.ExtFunIO.jointi.arm_r,model_info.ExtFunIO.jointi.arm_l];
orderArm = orderArm-min(orderArm)+1;
orderArmInv = [model_info.ExtFunIO.jointi.arm_l,model_info.ExtFunIO.jointi.arm_r];
orderArmInv = orderArmInv-min(orderArmInv)+1;

model_info.ExtFunIO.symQs.orderQs = orderQs;
model_info.ExtFunIO.symQs.orderQsInv = orderQsInv;
model_info.ExtFunIO.symQs.orderQsOpp1 = orderQsOpp1;
model_info.ExtFunIO.symQs.orderArm = orderArm;
model_info.ExtFunIO.symQs.orderArmInv = orderArmInv;

if strcmp(S.ModelName,'Rajagopal')
    % indexes to select kinematics left and right leg
    IndexLeft = model_info.ExtFunIO.jointi.leg_l;
    IndexRight = model_info.ExtFunIO.jointi.leg_r;
elseif strcmp(S.ModelName,'Gait92')
    % indexes to select kinematics left and right leg
    IndexLeft = [model_info.ExtFunIO.jointi.leg_l model_info.ExtFunIO.jointi.torso];
    IndexRight = [model_info.ExtFunIO.jointi.leg_r model_info.ExtFunIO.jointi.torso];
end

% indexes for symmetry steps
QsInvA = [model_info.ExtFunIO.coordi.pelvis_tilt,...
    model_info.ExtFunIO.coordi.pelvis_ty,...
    model_info.ExtFunIO.jointi.leg_r,model_info.ExtFunIO.jointi.leg_l,...
    model_info.ExtFunIO.coordi.lumbar_extension,...
    model_info.ExtFunIO.jointi.arm_r,model_info.ExtFunIO.jointi.arm_l]';
QsInvB = [model_info.ExtFunIO.coordi.pelvis_tilt,...
    model_info.ExtFunIO.coordi.pelvis_ty,...
    model_info.ExtFunIO.jointi.leg_l,model_info.ExtFunIO.jointi.leg_r,...
    model_info.ExtFunIO.coordi.lumbar_extension,...
    model_info.ExtFunIO.jointi.arm_l,model_info.ExtFunIO.jointi.arm_r]';

QdotsInvA = [model_info.ExtFunIO.coordi.pelvis_tilt,...
    model_info.ExtFunIO.coordi.pelvis_tx,model_info.ExtFunIO.coordi.pelvis_ty,...
    model_info.ExtFunIO.jointi.leg_r,model_info.ExtFunIO.jointi.leg_l,...
    model_info.ExtFunIO.coordi.lumbar_extension,...
    model_info.ExtFunIO.jointi.arm_r,model_info.ExtFunIO.jointi.arm_l]';
QdotsInvB = [model_info.ExtFunIO.coordi.pelvis_tilt,...
    model_info.ExtFunIO.coordi.pelvis_tx,model_info.ExtFunIO.coordi.pelvis_ty,...
    model_info.ExtFunIO.jointi.leg_l,model_info.ExtFunIO.jointi.leg_r,...
    model_info.ExtFunIO.coordi.lumbar_extension,...
    model_info.ExtFunIO.jointi.arm_l,model_info.ExtFunIO.jointi.arm_r]';

orderQsOpp = [model_info.ExtFunIO.coordi.pelvis_list:model_info.ExtFunIO.coordi.pelvis_list,...
    model_info.ExtFunIO.coordi.pelvis_rotation:model_info.ExtFunIO.coordi.pelvis_rotation,...
    model_info.ExtFunIO.coordi.pelvis_tz:model_info.ExtFunIO.coordi.pelvis_tz,...
    model_info.ExtFunIO.coordi.lumbar_bending:model_info.ExtFunIO.coordi.lumbar_bending,...
    model_info.ExtFunIO.coordi.lumbar_rotation:model_info.ExtFunIO.coordi.lumbar_rotation];

model_info.ExtFunIO.symQs.IndexLeft = IndexLeft;
model_info.ExtFunIO.symQs.IndexRight = IndexRight;
model_info.ExtFunIO.symQs.QsInvA = QsInvA;
model_info.ExtFunIO.symQs.QsInvB = QsInvB;
model_info.ExtFunIO.symQs.QdotsInvA = QdotsInvA;
model_info.ExtFunIO.symQs.QdotsInvB = QdotsInvB;
model_info.ExtFunIO.symQs.orderQsOpp = orderQsOpp;

model_info.muscle_info.IndexCalf = [model_info.muscle_info.muscle_index.lat_gas_r,...
    model_info.muscle_info.muscle_index.med_gas_r,...
    model_info.muscle_info.muscle_index.soleus_r,...
    model_info.muscle_info.muscle_index.lat_gas_l,...
    model_info.muscle_info.muscle_index.med_gas_l,...
    model_info.muscle_info.muscle_index.soleus_l];
end

