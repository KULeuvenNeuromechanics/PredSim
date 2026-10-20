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
% --------------------------------------------------------------------------
% This file is part of PredSim.
% 
% PredSim: A Framework for Rapid Predictive Simulations of Locomotion
% Copyright (c) 2026 KU Leuven
% 
% PredSim is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Affero General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% 
% PredSim is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public 
% License for more details.
% 
% You should have received a copy of the GNU Affero General Public License 
% along with PredSim. If not, see <https://www.gnu.org/licenses/>.
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

