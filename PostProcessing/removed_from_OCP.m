
%% removed from OCP_formulation, to be moved to post-processing

% Get muscle excitations from time derivative of muscle activations
e_opt_unsc = computeExcitationRaasch(a_opt_unsc,vA_opt_unsc,...
    ones(1,NMuscle)*tdeact,ones(1,NMuscle)*tact);

% Muscle excitations
vA_GC = zeros(N*2,NMuscle);
vA_GC(1:N-IC1i_c+1,:) = vA_opt_unsc(IC1i_c:end,:);
vA_GC(N-IC1i_c+2:N-IC1i_c+1+N,:) = vA_opt_unsc(1:end,orderMusInv);
vA_GC(N-IC1i_c+2+N:2*N,:) = vA_opt_unsc(1:IC1i_c-1,:);
% If the first heel strike was on the left foot then we invert so that
% we always start with the right foot, for analysis purpose
if strcmp(HS1,'l')
    vA_GC(:,:) = vA_GC(:,orderMusInv);
end
e_GC = computeExcitationRaasch(Acts_GC,vA_GC,...
    ones(1,NMuscle)*tdeact,ones(1,NMuscle)*tact);




JointAngle.labels = [{'time'}, fields(model_info.ExtFunIO.coordi)', {'pro_sup_l'},{'pro_sup_r'}];
% Two gait cycles
% Joint angles
q_opt_GUI_GC_2 = [q_opt_GUI_GC;q_opt_GUI_GC];
q_opt_GUI_GC_2(2*N+1:4*N,1) = q_opt_GUI_GC_2(2*N+1:4*N,1) + ...
    q_opt_GUI_GC_2(end,1) + ...
    q_opt_GUI_GC_2(end,1)-q_opt_GUI_GC_2(end-1,1);
q_opt_GUI_GC_2(2*N+1:4*N,model_info.ExtFunIO.coordi.pelvis_tx+1) = ...
    q_opt_GUI_GC_2(2*N+1:4*N,model_info.ExtFunIO.coordi.pelvis_tx+1) + ...
    2*q_opt_unsc_all.deg(end,model_info.ExtFunIO.coordi.pelvis_tx);
% Muscle activations (to have muscles turning red when activated).
Acts_GC_GUI = [Acts_GC;Acts_GC];
% Combine data joint angles and muscle activations
JointAngleMuscleAct.data = [q_opt_GUI_GC_2,Acts_GC_GUI];
% Get muscle labels
muscleNamesAll = model_info.muscle_info.muscle_names;
% Combine labels joint angles and muscle activations
JointAngleMuscleAct.labels = JointAngle.labels;
for i = 1:NMuscle
    JointAngleMuscleAct.labels{i+size(q_opt_GUI_GC_2,2)} = ...
        [muscleNamesAll{i},'/activation'];
end
JointAngleMuscleAct.inDeg = 'yes';
filenameJointAngles = fullfile(OutFolder,[S.subject.name '.mot']);
write_motionFile_v40(JointAngleMuscleAct, filenameJointAngles);