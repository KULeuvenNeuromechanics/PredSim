function [] = PostProcessing(S,model_info,f_casadi)
% --------------------------------------------------------------------------
% PostProcessing
%   This function calls subfunctions that post-process the simulation
%   results.
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
%   - f_casadi -
%   * Struct containing all casadi functions.
% 
% Original author: Lars D'Hondt
% Original date: May/2022
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

% load results
Outname = fullfile(S.misc.save_folder,[S.misc.result_filename '.mat']);
load(Outname,'R');

% kinematics in radians
R.kinematics.Qs_rad = R.kinematics.Qs;
R.kinematics.Qs_rad(:,model_info.ExtFunIO.jointi.rotations) = ...
    R.kinematics.Qs(:,model_info.ExtFunIO.jointi.rotations)*pi/180;
R.kinematics.Qdots_rad = R.kinematics.Qdots;
R.kinematics.Qdots_rad(:,model_info.ExtFunIO.jointi.rotations) = ...
    R.kinematics.Qdots(:,model_info.ExtFunIO.jointi.rotations)*pi/180;
R.kinematics.Qddots_rad = R.kinematics.Qddots;
R.kinematics.Qddots_rad(:,model_info.ExtFunIO.jointi.rotations) = ...
    R.kinematics.Qddots(:,model_info.ExtFunIO.jointi.rotations)*pi/180;

% misc info for quick access
R.misc.body_mass = model_info.mass;
R.misc.body_weight = model_info.mass*9.81;

%% Call post-processing subfunctions
% Do mind the order of the functions, since some of them use information
% computed by another function.

% Saves a .mot file containing 2 full gait cycles, for visualisation in the
% OpenSim GUI.
[R] = PostProcess_write_motion_file(model_info,f_casadi,R);

% Get Inverse Dynamic torque (or force) for each coordinate.
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

% Evaluate ligament lengths, forces, powers
[R] = PostProcessing_ligaments(model_info,f_casadi,R);

% Calculate all terms of the passive coordinate moments.
[R] = PostProcess_passive_moments(model_info,f_casadi,R);

% Compute spatio-temporal characteristics of the full gait cycle.
[R] = PostProcess_spatio_temporal(model_info,f_casadi,R);

% Evaluate the metabolic energetics for all implemented metabolic energy
% models.
[R] = PostProcess_metabolic_energy(model_info,f_casadi,R);

% Calculate orthosis forces
[R] = PostProcess_orthosis(model_info,f_casadi,R);

% Please implement additional post-processing steps as functions following
% the template, and call them from here.
% [R] = PostProcessing_subfunction_template(model_info,f_casadi,R);


%%
save(Outname,'R','-append');


