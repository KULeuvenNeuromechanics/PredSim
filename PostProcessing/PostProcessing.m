function [] = PostProcessing(S,model_info,f_casadi)


%% load results
Outname = fullfile(S.subject.save_folder,[S.post_process.result_filename '.mat']);
load(Outname,'R');

%%

t_mesh = [R.t_mesh(1:end-1),R.t_mesh(1:end-1)+R.t_mesh(end),...
    R.t_mesh(1:end-1)+R.t_mesh(end)*2,R.t_mesh(1:end-1)+R.t_mesh(end)*3];
%%


%% Create .mot file for OpenSim GUI
writeIKmotion = 1;
if writeIKmotion
    q_opt_GUI_GC_1 = [R.Qs];
    q_opt_GUI_GC_1(:,model_info.ExtFunIO.jointi.rotations) =...
        q_opt_GUI_GC_1(:,model_info.ExtFunIO.jointi.rotations)*pi/180;
    q_opt_GUI_GC_2 = q_opt_GUI_GC_1;
    q_opt_GUI_GC_2(:,model_info.ExtFunIO.jointi.base_forward) =...
        q_opt_GUI_GC_2(:,model_info.ExtFunIO.jointi.base_forward) +...
        q_opt_GUI_GC_1(end,model_info.ExtFunIO.jointi.base_forward);
    JointAngle.labels = [{'time'},model_info.ExtFunIO.coord_names.all(:)'];
    % Two gait cycles
    % Joint angles
    q_opt_GUI_GC = [t_mesh',[q_opt_GUI_GC_1;q_opt_GUI_GC_2]];
    % Muscle activations (to have muscles turning red when activated).
    Acts_GC = R.a;
    Acts_GC_GUI = [Acts_GC;Acts_GC];
    % Combine data joint angles and muscle activations
    JointAngleMuscleAct.data = [q_opt_GUI_GC,Acts_GC_GUI];
    % Combine labels joint angles and muscle activations
    JointAngleMuscleAct.labels = JointAngle.labels;
    for i = 1:model_info.muscle_info.NMuscle
        JointAngleMuscleAct.labels{i+size(q_opt_GUI_GC,2)} = ...
            [model_info.muscle_info.muscle_names{i},'/activation'];
    end
    JointAngleMuscleAct.inDeg = 'true';
%     OutFolder = fullfile(pathRepo,'Results',S.ResultsFolder);
    filenameJointAngles = fullfile(S.subject.save_folder,[S.post_process.result_filename '.mot']);
    write_motionFile_v40(JointAngleMuscleAct, filenameJointAngles);
end

