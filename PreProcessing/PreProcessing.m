function [S,model_info] = PreProcessing(S,osim_path)
%
% Prepare the information that is required for the problem formulation.
%
% Author: Lars D'Hondt
%
% Date: 18/January/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% subject specific folder
% to store information about model
if ~isfolder(S.subject.save_folder)
    mkdir(S.subject.save_folder);
end

%% generate external function
% skeletal and contact dynamics, compiled based on .osim file
% https://github.com/Lars-DHondt-KUL/opensimAD
ext_fun_path = fullfile(S.misc.subject_folder,['F_' osim_file_name '.dll']);
if ~exist(ext_fun_path,'file')
    osim2dll(S,osim_path);
end

%% create model_info
% book-keeping
model_info = get_model_info(S,osim_path);

%% read Muscle-Tendon parameters, and scale them if needed
% lMo, vMo, FMo, alphao, lTs
model_info = read_and_scale_MTparameters(S,osim_path,model_info);

%% analyse musculo-skeletal geometry, and fit an expression
% based on OpenSim muscle analysis tool
get_musculoskeletal_geometry_approximation(S,osim_path,model_info);

%% create additional fields with information that are needed later
%
model_info = update_model_info(model_info);


