function [] = copy_modelsettingsfile(or_model_path,output_model_path)
%copy_musclegeom_information copies muscle geometry information file from
%or_model_path to output_model_path
%   Detailed explanation goes here

% folder original model
[folder_or,modelname,~] = fileparts(or_model_path);
settinsfile_or = fullfile(folder_or,['settings_' modelname '.m']);
if exist(settinsfile_or,'file')

    % folder output
    [folder_output,modelname,~] = fileparts(output_model_path);
    settingsfile_out = fullfile(folder_output,['settings_' modelname '.m']);

    % copy file
    copyfile(settinsfile_or,settingsfile_out);
else
    disp(['No settings file for model :', modelname]);
end





end