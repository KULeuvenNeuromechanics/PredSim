function [S,model_info] = PreProcessing(S,osim_path)
% --------------------------------------------------------------------------
% PreProcessing
%   This function calls all preprocessing steps needed to prepare formulating
%   the OCP.
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%
% OUTPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% Original author: Lars D'Hondt
% Original date: 7/March/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


% Create external function to describe the rigid-body skeletal dynamics and foot-ground contact dynamics.
osim2dll(S,osim_path);

% Create a struct to contain all information about the neuro-musculoskeletal model
model_info = get_model_info(S,osim_path);

% Finalize model_info
model_info = update_model_info(S,model_info);

% Read muscle-tendon parameters from the opensim model, and scale them 
model_info = read_and_scale_MTparameters(S,osim_path,model_info);

% Read actuator parameters from opensim model
model_info = get_actuator_info(S,osim_path,model_info);

% Describe uscle-tendon lengths, velocities, and moment arms in function of coordinate values
model_info = get_musculoskeletal_geometry_approximation(S,osim_path,model_info);




