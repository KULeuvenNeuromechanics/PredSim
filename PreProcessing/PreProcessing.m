function [S,model_info] = preprocessing(S,osim_path)
%
% This function calls all preprocessing steps needed to prepare formulating
% the OCP.
%
% Author: Lars D'Hondt
% Date: 7/March/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create external function to describe the rigid-body skeletal dynamics and foot-ground contact dynamics.
osim2dll(S,osim_path);

% Create a struct to contain all information about the neuro-musculoskeletal model
model_info = get_model_info(S,osim_path);

% Read muscle-tendon parameters from the opensim model, and scale them 
model_info = read_and_scale_MTparameters(S,osim_path,model_info);




model_info = get_musculoskeletal_geometry_approximation(S,osim_path,model_info);

model_info = update_model_info(model_info);
