function [] = completeModelScaling(osim_scaleTool, osim_generic, opts)
% --------------------------------------------------------------------------
% completeModelScaling
%   After running the OpenSim scale tool, run this function to scale model
%   parameters that were not yet scaled but are relevant for PredSim.
%   Scale factors are calculated from body height and mass, based on [1].
%   Contact spheres are scaled as in [2].
%   
%   | Parameter                  | Scaling       | Unscaled values from |
%   |----------------------------|---------------|----------------------|
%   | muscle max isometric force | mass^(2/3)    | osim_scaleTool       |
%   | ligament max force         | mass^(2/3)    | osim_scaleTool       |
%   | actuator max torque        | mass*height   | osim_generic         |
%   | contact sphere radius      | foot size     | osim_scaleTool       |
%   | contact sphere position    | foot size     | osim_scaleTool       |
%   | contact sphere stiffness   | mass/height^2 | osim_scaleTool       |
%   | contact sphere dissipation | 1/height      | osim_scaleTool       |
% 
%
%   This function also saves a .mat file in the same folder as osim_scaled
%   (osim_filename + '_scaling'). This file contains a struct with scale
%   factors for PredSim settings that are not size-independent:
%   - length            distance constraints
%   - moment            stiffness and damping coefficients of rotational DOFs
%   - weight_Emetab     weight factor of metabolic energy objective term [2]
%   - weight_pass_torq  weight factor of passive torque objective term [2]
%
%
%   References
%   [1] A. L. Hof, “Scaling gait data to body size,” Gait & Posture, 
%   vol. 4, no. 3, pp. 222–223, May 1996, doi: 10.1016/0966-6362(95)01057-2.
%   [2] I. Vandekerckhove et al., “Muscle weakness but also contractures 
%   contribute to the progressive gait pathology in children with Duchenne 
%   muscular dystrophy: a simulation study,” J NeuroEngineering Rehabil, 
%   vol. 22, no. 1, p. 103, May 2025, doi: 10.1186/s12984-025-01631-x.
%
%
%   See also
%   AdaptOpenSimModel
%   scaleMuscleForce scaleLigaments 
%   scaleContactSpheres 
%   scaleActuators
%
%
% INPUT:
%
% Required:
%   - osim_scaleTool -
%   * Path to the OpenSim model file created by the OpenSim Scale Tool.
%   (i.e. input model file)
% 
%   - osim_generic -
%   * Path to the OpenSim model file with the unscaled model. This model
%   serves as reference for actuators and mass.
%
% Optional:
% (name-value pairs)
%
%   - osim_scaled -  Default: osim_scaleTool
%   * Path to the OpenSim model file with the fully scaled model. (i.e.
%   output model file)
%
%   - height_subject -  Default: 1.65
%   * Height of the subject (i.e. scaled model) in meter. Used to
%   calculate height ratio for scaling.
%
%   - height_generic -  Default: 1.65
%   * Height of the generic (unscaled) model in meter. Used to
%   calculate height ratio for scaling. Default value of 1.65m is for
%   models used in Falisse et al. (2022) and D'Hondt et al. (2024).
%
%   - scale_factor_foot_left - Default: 1
%   * Scale factor that was used to scale left foot in OpenSim scale tool.
%
%   - scale_factor_foot_right - Default: scale_factor_foot_left
%   * Scale factor that was used to scale right foot in OpenSim scale tool.
%
% Example:
% completeModelScaling(osim_scaleTool, osim_generic,...
%     height_subject=1.72, ...
%     scale_factor_foot_left=0.5, ...
%     scale_factor_foot_right=0.5)
%
%   The model name is updated with the following suffixes:
%   - _AdaptedForPredSim     indicates the model has been adapted for PredSim
%   - _sf                    indicates muscle forces have been scaled.
%   A warning is shown if either suffix is already present.
% 
% Original author: Lars D'Hondt
% Original date: 3 January 2025

% Last edit by: Ellis Van Can
%  Last edit date: 5 May 2026
% --------------------------------------------------------------------------
arguments
    osim_scaleTool (1,:) char
    osim_generic   (1,:) char
    opts.osim_scaled             (1,:) char   = '';
    opts.height_subject          (1,1) double = 1.65;
    opts.height_generic          (1,1) double = 1.65;
    opts.scale_factor_foot_left  (1,1) double = 1;
    opts.scale_factor_foot_right (1,1) double = nan;  
end

% handle defaults that depend on other arguments
if isempty(opts.osim_scaled)
    opts.osim_scaled = osim_scaleTool;
end
if isnan(opts.scale_factor_foot_right)
    opts.scale_factor_foot_right = opts.scale_factor_foot_left;
end

% unpack opts
osim_scaled             = opts.osim_scaled;
height_subject          = opts.height_subject;
height_generic          = opts.height_generic;
scale_factor_foot_left  = opts.scale_factor_foot_left;
scale_factor_foot_right = opts.scale_factor_foot_right;

% path to helper functions
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath(fullfile(pathRepo, 'VariousFunctions'))

%%
mass_generic = getModelMass(osim_generic);
mass_subject = getModelMass(osim_scaleTool);

sf_mass = mass_subject/mass_generic;
sf_length = height_subject/height_generic;

sf_moment = sf_mass * sf_length;
sf_contact_stiffness = 1; %sf_mass/sf_length^2;
sf_contact_dissipation = 1; %1/sf_length;

scale_factors.mass = sf_mass;
scale_factors.length = sf_length;

scale_factors.moment = sf_moment;

scale_factors.weight_Emetab = 1/(sf_mass*sf_length^2);
scale_factors.weight_pass_torq = 1/sf_moment^2;

% save scale factors for later use i.e. to scale parameters and cost
% function weights via PredSim Settings
save(replace(osim_scaled,'.osim','_scaling.mat'), "scale_factors");

%%

scaleMuscleForce(osim_scaled, mass_generic, osim_scaleTool)

scaleLigaments(osim_scaled, mass_generic)

%% coordinate actuators

scaleActuators(osim_scaled, sf_moment)


%% contact spheres
sf_contact.stiffness = sf_contact_stiffness;
sf_contact.dissipation = sf_contact_dissipation;
sf_contact.foot_left = scale_factor_foot_left;
sf_contact.foot_right = scale_factor_foot_right;

scaleContactSpheres(osim_scaled, osim_scaled, sf_contact)

%% Check scaling history and update model name
import org.opensim.modeling.*;
model = Model(osim_scaleTool);
model_name = char(model.getName);

if contains(model_name, 'AdaptedForPredSim')
    warning('Model "%s" has already been adapted for PredSim.', model_name);
else
    model.setName([model_name '_AdaptedForPredSim']);
end

if contains(model_name, '_sf')                        
    warning('Model "%s" has already been force-scaled.', model_name);
else
    model.setName([char(model.getName) '_sf']);         
end


model.print(osim_scaled);
fprintf("The scaled OpenSim model is saved as:\n\t'%s'\n", osim_scaled)

end % end of function
