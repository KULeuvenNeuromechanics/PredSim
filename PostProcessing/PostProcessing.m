function [] = PostProcessing(S,model_info,f_casadi)


%% load results
Outname = fullfile(S.subject.save_folder,[S.post_process.result_filename '.mat']);
load(Outname,'R');

%% Create .mot file for OpenSim GUI
[R] = PostProcess_write_motion_file(S,model_info,f_casadi,R);


%% Add kinematics in radians
R.Qs_rad(:,model_info.ExtFunIO.jointi.rotations) = R.Qs(:,model_info.ExtFunIO.jointi.rotations)*pi/180;
R.Qdots_rad(:,model_info.ExtFunIO.jointi.rotations) = R.Qdots(:,model_info.ExtFunIO.jointi.rotations)*pi/180;
R.Qddots_rad(:,model_info.ExtFunIO.jointi.rotations) = R.Qddots(:,model_info.ExtFunIO.jointi.rotations)*pi/180;

%% 
[R] = PostProcess_external_function(S,model_info,f_casadi,R);




