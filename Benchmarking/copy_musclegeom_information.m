function [] = copy_musclegeom_information(or_model_path,output_model_path)
%copy_musclegeom_information copies muscle geometry information file from
%or_model_path to output_model_path
%   Detailed explanation goes here

% folder original model
[folder_or,modelname,~] = fileparts(or_model_path);
geom_file = fullfile(folder_or,[modelname '_f_lMT_vMT_dM_poly_3_9']);

% folder output
[folder_output,modelname,~] = fileparts(output_model_path);
geom_file_out = fullfile(folder_output,[modelname '_f_lMT_vMT_dM_poly_3_9']);


% copy file
copyfile(geom_file,geom_file_out);



end