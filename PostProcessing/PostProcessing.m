function [] = PostProcessing(S,model_info,f_casadi)


% load results
Outname = fullfile(S.subject.save_folder,[S.post_process.result_filename '.mat']);
load(Outname,'R');

% kinematics in radians
R.kinematics.Qs_rad = R.kinematics.Qs;
R.kinematics.Qs_rad(:,model_info.ExtFunIO.jointi.rotations) = R.kinematics.Qs(:,model_info.ExtFunIO.jointi.rotations)*pi/180;
R.kinematics.Qdots_rad = R.kinematics.Qdots;
R.kinematics.Qdots_rad(:,model_info.ExtFunIO.jointi.rotations) = R.kinematics.Qdots(:,model_info.ExtFunIO.jointi.rotations)*pi/180;
R.kinematics.Qddots_rad = R.kinematics.Qddots;
R.kinematics.Qddots_rad(:,model_info.ExtFunIO.jointi.rotations) = R.kinematics.Qddots(:,model_info.ExtFunIO.jointi.rotations)*pi/180;

% misc info for quick access
R.misc.body_mass = S.subject.mass;
R.misc.body_weight = S.subject.mass*9.81;

%% Call post-processing subfunctions
% Do mind the order of the functions, since some of them use information
% computed by another function.

% Saves a .mot file containing 2 full gait cycles, for visualisation in the
% OpenSim GUI.
[R] = PostProcess_write_motion_file(model_info,f_casadi,R);

% Get Inverse Dynamic torque (or force) for each cordinate.
[R] = PostProcess_get_ID(model_info,f_casadi,R);

% Compute variables related to foot-ground contact.
[R] = PostProcess_ground_reaction(model_info,f_casadi,R);

% Reconstruct muscle excitations from implicit activation dynamics
% formulation.
[R] = PostProcess_muscle_excitations(model_info,f_casadi,R);

% Evaluate the approximated muscle-tendon lenghts and -velocities, and
% moment arms. Compare the approximated with the geometry from the original
% OpenSim model.
[R] = PostProcess_msk_geometry(model_info,f_casadi,R);

% Reconstruct all force- and length components of muscle fiber and tendon.
[R] = PostProcess_muscletendon_dynamics(model_info,f_casadi,R);

% Calculate all terms of the passive coordinate moments.
[R] = PostProcess_passive_moments(model_info,f_casadi,R);

% Compute spatio-temporal characteristics of the full gait cycle.
[R] = PostProcess_spatio_temporal(model_info,f_casadi,R);

% Evaluate the metobiloc energetics for all implemented metabolic energy
% models.
[R] = PostProcess_metabolic_energy(model_info,f_casadi,R);

% Please implement additional post-processing steps as functions following
% the template, and call them from here.
% [R] = PostProcessing_subfunction_template(model_info,f_casadi,R);


%%
Outname = fullfile(S.subject.save_folder,[S.post_process.result_filename '.mat']);
load(Outname,'w_opt','stats','setup','model_info');
save(Outname,'w_opt','stats','setup','R','model_info');


