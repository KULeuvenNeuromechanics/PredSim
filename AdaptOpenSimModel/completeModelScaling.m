function [] = completeModelScaling(osim_scaleTool, osim_generic, osim_scaled,...
    height_subject, height_generic, scale_factor_foot_left, scale_factor_foot_right)
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
%   | contact sphere radius      | foot          | osim_generic         |
%   | contact sphere position    | foot          | osim_generic         |
%   | contact sphere stiffness   | mass*height^2 | osim_generic         |
%   | contact sphere dissipation | height        | osim_generic         |
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
%   [2] I. Vandekerckhove et al (unpublished)
%
%   See also 
%   AdaptOpenSimModel
%   scaleMuscleForce scaleLigaments 
%   get_contact_spheres removeContactSpheres add_contact_spheres 
%   getActuators removeActuators add_actuators
%
%
% INPUT:
%   - osim_scaleTool -
%   * Path to the OpenSim model file created by the OpenSim Scale Tool.
% 
%   - osim_generic -
%   * Path to the OpenSim model file with the unscaled model. This model
%   serves as reference for contact spheres and actuators.
%
%   - osim_scaled - (optional) Default: osim_scaleTool
%   * Path to the OpenSim model file with the fully scaled model.
%
%   - height_subject - (optional) Default: 1.65
%   * Height of the subject (i.e. scaled model) in meter. Used to
%   calculated height ratio for scaling.
%
%   - height_generic - (optional) Default: 1.65
%   * Height of the generic (unscaled) model in meter. Used to
%   calculate height ratio for scaling. Default value of 1.65m is for
%   models used in Falisse et al. (2022) and D'Hondt et al. (2024).
%
%   - scale_factor_foot_left - (optional) Default: 1
%   * Scale factor that was used to scale left foot in OpenSim  scale tool.
%
%   - scale_factor_foot_right - (optional) Default: scale_factor_foot_left
%   * Scale factor that was used to scale right foot in OpenSim  scale tool.
%
% 
% Original author: Lars D'Hondt
% Original date: 3 January 2025
% --------------------------------------------------------------------------

arguments
    osim_scaleTool
    osim_generic
    osim_scaled = osim_scaleTool;
    height_subject = 1.65;
    height_generic = 1.65;
    scale_factor_foot_left = 1;
    scale_factor_foot_right = scale_factor_foot_left;
end

addpath('../VariousFunctions')

%%

mass_generic = getModelMass(osim_generic);
mass_subject = getModelMass(osim_scaleTool);

sf_mass = mass_subject/mass_generic;
sf_length = height_subject/height_generic;

sf_moment = sf_mass * sf_length;
sf_contact_stiffness = sf_mass*sf_length^2;
sf_contact_dissipation = sf_length;

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

% get generic values
[torq_act] = getActuators(osim_generic);

% scale values
for i=1:length(torq_act)
    torq_act(i).max_torque = torq_act(i).max_torque * sf_moment;
end

% remove actuators from scaled model before adding
removeActuators(osim_scaled)

% add to scaled model
add_actuators(osim_scaled,torq_act)

%% contact spheres

% get generic values
[contact_spheres] = get_contact_spheres(osim_generic);

% scale values
for i=1:length(contact_spheres)

    contact_spheres(i).stiffness = contact_spheres(i).stiffness*sf_contact_stiffness;
    contact_spheres(i).dissipation = contact_spheres(i).dissipation*sf_contact_dissipation;

    name = contact_spheres(i).name;
    [leftname,~] = mirrorName(name);
    if strcmp(name,leftname)
        contact_spheres(i).radius = contact_spheres(i).radius*scale_factor_foot_left(1);
        contact_spheres(i).location = contact_spheres(i).location.*scale_factor_foot_left;
    else
        contact_spheres(i).radius = contact_spheres(i).radius*scale_factor_foot_right(1);
        contact_spheres(i).location = contact_spheres(i).location.*scale_factor_foot_right;
    end
end

% remove contact spheres from scaled model to avoid duplicating them
removeContactSpheres(osim_scaled)

% add scaled contact spheres
add_contact_spheres(osim_scaled,contact_spheres)


end % end of function