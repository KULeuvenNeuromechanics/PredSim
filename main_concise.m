clear
close all
clc
% path to the repository folder
[pathRepo,~,~] = fileparts(mfilename('fullpath'));
% path to the folder that contains the repository folder
[pathRepoFolder,~,~] = fileparts(pathRepo);


addpath(fullfile(pathRepo,'DefaultSettings'))


model_name = 'DHondt_et_al_2024_3seg';


% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',model_name,[model_name '.osim']);


[S] = initializeSettings(model_name);

S.misc.main_path = pathRepo;

% path to folder where you want to store the results of the OCP
S.subject.save_folder  = fullfile(pathRepoFolder,'PredSimResults',S.subject.name);

S.solver.CasADi_path = casadi.GlobalOptions.getCasadiPath(); % ask casadi

S.OpenSimADOptions.verbose_mode = 0;

S.OpenSimADOptions.jointsOrder = {'ground_pelvis', 'hip_l', 'hip_r', 'knee_l',...
    'knee_r', 'ankle_l', 'ankle_r', 'subtalar_l', 'subtalar_r', 'midtarsal_l',...
    'midtarsal_r', 'tarsometatarsal_l', 'tarsometatarsal_r', 'mtp_l', 'mtp_r',...
    'back', 'acromial_l', 'acromial_r', 'elbow_l', 'elbow_r', 'radioulnar_l',...
    'radioulnar_r', 'radius_hand_l', 'radius_hand_r'};
S.OpenSimADOptions.coordinatesOrder = {
     'pelvis_tilt', 'pelvis_list', 'pelvis_rotation', 'pelvis_tx', 'pelvis_ty',... 
     'pelvis_tz', 'hip_flexion_l', 'hip_adduction_l', 'hip_rotation_l',... 
     'hip_flexion_r', 'hip_adduction_r', 'hip_rotation_r', 'knee_angle_l',... 
     'knee_angle_r', 'ankle_angle_l', 'ankle_angle_r', 'subtalar_angle_l',...
     'subtalar_angle_r', 'mtj_angle_l', 'mtj_angle_r', 'mtp_angle_l',...
     'mtp_angle_r', 'lumbar_extension', ...
     'lumbar_bending', 'lumbar_rotation', 'arm_flex_l', 'arm_add_l', ...
     'arm_rot_l', 'arm_flex_r', 'arm_add_r', 'arm_rot_r', 'elbow_flex_l', ...
     'elbow_flex_r'};

% S.solver.max_iter       = 5;

S.solver.N_meshes = 50;

S.solver.run_as_batch_job = 1;



% Start simulation
if S.solver.run_as_batch_job
    add_pred_sim_to_batch(S,osim_path)
else
    run_pred_sim(S,osim_path);
end

