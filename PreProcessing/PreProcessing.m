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
% --------------------------------------------------------------------------
% This file is part of PredSim.
% 
% PredSim: A Framework for Rapid Predictive Simulations of Locomotion
% Copyright (c) 2026 KU Leuven
% 
% PredSim is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Affero General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% 
% PredSim is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public 
% License for more details.
% 
% You should have received a copy of the GNU Affero General Public License 
% along with PredSim. If not, see <https://www.gnu.org/licenses/>.
% --------------------------------------------------------------------------

% Update settings about points that need to be exported from external function
[S] = updateExport3DPositionsVelocities(S,osim_path);

% Test the user-defined orthosis functions and add the required OpenSimAD
% options
[S] = finaliseOrthosisDefinitions(S,osim_path);

% Create external function to describe the rigid-body skeletal dynamics and
% foot-ground contact dynamics.
S = osim2dll(S,osim_path);

% Create a struct to contain all information about the neuro-musculoskeletal model
model_info = get_model_info(S,osim_path);

% Read ligament parameters from the opensim model and settings
if model_info.ligament_info.NLigament > 0
    model_info = get_ligament_info(S,osim_path,model_info);
end

% Read muscle-tendon parameters from the opensim model, and scale them 
model_info = read_and_scale_MTparameters(S,osim_path,model_info);

% Read actuator parameters from opensim model
model_info = get_actuator_info(S,osim_path,model_info);

% Describe muscle-tendon lengths, velocities, and moment arms in function of coordinate values
model_info = get_musculoskeletal_geometry_approximation(S,osim_path,model_info);

% Read actuator parameters from opensim model
model_info = get_passive_moment_info(S,model_info);

% Finalize model_info
model_info = update_model_info(S,osim_path,model_info);


end
