function [] = PostProcessing(S,model_info,f_casadi)


% load results
Outname = fullfile(S.subject.save_folder,[S.post_process.result_filename '.mat']);
load(Outname,'R');

% Add kinematics in radians
R.kinematics.Qs_rad(:,model_info.ExtFunIO.jointi.rotations) = R.kinematics.Qs(:,model_info.ExtFunIO.jointi.rotations)*pi/180;
R.kinematics.Qdots_rad(:,model_info.ExtFunIO.jointi.rotations) = R.kinematics.Qdots(:,model_info.ExtFunIO.jointi.rotations)*pi/180;
R.kinematics.Qddots_rad(:,model_info.ExtFunIO.jointi.rotations) = R.kinematics.Qddots(:,model_info.ExtFunIO.jointi.rotations)*pi/180;

%%
[R] = PostProcess_write_motion_file(S,model_info,f_casadi,R);

[R] = PostProcess_get_ID(S,model_info,f_casadi,R);

[R] = PostProcess_ground_reaction(S,model_info,f_casadi,R);

[R] = PostProcess_muscle_excitations(S,model_info,f_casadi,R);

[R] = PostProcess_msk_geometry(S,model_info,f_casadi,R);

[R] = PostProcess_muscletendon_dynamics(S,model_info,f_casadi,R);

[R] = PostProcess_metabolic_energy(S,model_info,f_casadi,R);

[R] = PostProcess_passive_moments(S,model_info,f_casadi,R);

[R] = PostProcess_spatio_temporal(S,model_info,f_casadi,R);

%%
Outname = fullfile(S.subject.save_folder,[S.post_process.result_filename '.mat']);
load(Outname,'w_opt','stats','setup','model_info');
save(Outname,'w_opt','stats','setup','R','model_info');


