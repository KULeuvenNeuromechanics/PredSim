function [R] = PostProcess_write_motion_file(model_info,f_casadi,R)
% --------------------------------------------------------------------------
% PostProcess_write_motion_file
%   This function creates a motionfile with 2 steps of predicted gait.
% 
% INPUT:
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
% Original date: 12/May/2022
%
% Last edit by: Bram Van Den Bosch
% Last edit date: 23/Sept/2024
% --------------------------------------------------------------------------

% Two gait cycles
t_mesh = [R.time.mesh_GC(1:end-1),R.time.mesh_GC(1:end-1)+R.time.mesh_GC(end)];
% Joint angles
q_opt_GUI_GC_1 = [R.kinematics.Qs];
q_opt_GUI_GC_2 = q_opt_GUI_GC_1;
q_opt_GUI_GC_2(:,model_info.ExtFunIO.jointi.base_forward) =...
    q_opt_GUI_GC_2(:,model_info.ExtFunIO.jointi.base_forward) +...
    R.spatiotemp.dist_trav;
JointAngle.labels = [{'time'},model_info.ExtFunIO.coord_names.all(:)'];

q_opt_GUI_GC = [t_mesh',[q_opt_GUI_GC_1;q_opt_GUI_GC_2]];
% Muscle activations (to have muscles turning red when activated).
Acts_GC = R.muscles.a;
Acts_GC_GUI = [Acts_GC;Acts_GC];
% Combine data joint angles and muscle activations
JointAngleMuscleAct.data = [q_opt_GUI_GC,Acts_GC_GUI];
% Combine labels joint angles and muscle activations
JointAngleMuscleAct.labels = JointAngle.labels;
for i = 1:model_info.muscle_info.NMuscle
    JointAngleMuscleAct.labels{i+size(q_opt_GUI_GC,2)} = ...
        [model_info.muscle_info.muscle_names{i},'/activation'];
end
JointAngleMuscleAct.inDeg = 'yes';
filenameJointAngles = fullfile(R.S.misc.save_folder,...
    [R.S.misc.result_filename '.mot']);
write_motionFile_v40(JointAngleMuscleAct, filenameJointAngles);

% if gravity vector was tilted (walking on a slope) also export new mot
% file with gravity vector along y-axis.
if isfield(model_info,'slope') && abs(model_info.slope)>0
    % convert slope to angle
    fi = atan2(model_info.slope,1);

    % rotate data mot file
    data = JointAngleMuscleAct.data;
	Rotm =    [cos(fi)	-sin(fi)	0	
		sin(fi)	cos(fi)	0	
		0	0	1];

    i_pelvis_tilt = strcmp(JointAngleMuscleAct.labels,'pelvis_tilt');
    data(:,i_pelvis_tilt) = data(:,i_pelvis_tilt)+fi*180/pi;

    i_pelvis_tx = find(strcmp(JointAngleMuscleAct.labels,'pelvis_tx'));
    data_tpelvis =data(:,i_pelvis_tx:i_pelvis_tx+2);

    data(:,i_pelvis_tx) = data_tpelvis*Rotm(:,1);
    data(:,i_pelvis_tx+1) = data_tpelvis*Rotm(:,2);
    data(:,i_pelvis_tx+2) = data_tpelvis*Rotm(:,3);

    % structure for output
    slope_mot.data = data;
    slope_mot.labels = JointAngleMuscleAct.labels;
    slope_mot.inDeg = 'yes';

    % export mot file
    filenameJointAngles = fullfile(R.S.misc.save_folder,...
    [R.S.misc.result_filename '_slope.mot']);
    write_motionFile_v40(slope_mot, filenameJointAngles);


end

