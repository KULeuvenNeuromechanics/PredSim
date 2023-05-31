function [model_info] = read_and_scale_MTparameters(S,osim_path,model_info)
% --------------------------------------------------------------------------
% read_and_scale_MTparameters
%   Read muscle-tendon parameter values from the specified osim file:
%       - FMo: maximum isometric force
%       - lMo: optimal fiber length
%       - lTs: tendon slack length
%       - alphao: pennation angle at optimal fiber length
%       - vMmax: maximum contraction velocity
%       The obtained values can be scaled by a factor given in
%       "S.subject.MT_params"
%
%   Parameter values not obtained from the osim file:
%       - aTendon: tendon stiffness of the generic stress-strain relation.
%       Default value = 35. This value can be scaled by a
%       factor given in "S.subject.tendon_stiff"
%       - tensions: specific muscle tension. Data from Uchida et al. (2016)
%       is loaded through getSpecificTensions.m
%       - pctsts: fraction of muscle fibers that are slow twitch fibers. Data
%       from Uchida et al. (2016) is loaded through getSlowTwitchRatios.m
%       - muscle_strength: scale factor on the active muscle fiber force. The
%       default value of 1 can be replaced through "S.subject.muscle_strength"
%       - muscle_stiffness: positioning of passive force-length curve of a
%       muscle w.r.t. its normalized fiber length. Default value is 1, can
%       be changed to the value specified in "S.subject.muscle_stiff".
%
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% Original author: Lars D'Hondt
% Original date: 18/March/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

muscle_info = model_info.muscle_info;
muscleNames = muscle_info.muscle_names;
NMuscle = muscle_info.NMuscle;

%% read default muscle- and tendon parameters from opensim model file
[FMo, lMo, lTs, alphao, vMmax, fiber_damping, muscle_pass_stiff_scale,...
    tendon_stiff] = getMTparameters(osim_path,muscleNames);

%% parameters not from osim file
% strength: scale factor for active fiber force
muscle_strength = ones(1,NMuscle);

% stiffness shift: onset of passive fiber force-length w.r.t. norm fiber force
muscle_pass_stiff_shift = ones(1,NMuscle);


% specific tensions of muscle fibers (for metabolics)
specific_tension = getSpecificTensions(muscleNames)';

% ratio of slow twitch muscle fibers (for metabolics)
slow_twitch_fiber_ratio = getSlowTwitchRatios(muscleNames)';


%% Organise in struct

for i=1:NMuscle
    parameters(i).muscle_name = muscleNames{i};
end
[parameters] = double_array_to_struct_array(parameters,'FMo',FMo);
[parameters] = double_array_to_struct_array(parameters,'lMo',lMo);
[parameters] = double_array_to_struct_array(parameters,'lTs',lTs);
[parameters] = double_array_to_struct_array(parameters,'alphao',alphao);
[parameters] = double_array_to_struct_array(parameters,'vMmax',vMmax);
[parameters] = double_array_to_struct_array(parameters,'tendon_stiff',tendon_stiff);
[parameters] = double_array_to_struct_array(parameters,'tendon_stiff_shift',[]);
[parameters] = double_array_to_struct_array(parameters,'specific_tension',specific_tension);
[parameters] = double_array_to_struct_array(parameters,'slow_twitch_fiber_ratio',slow_twitch_fiber_ratio);
[parameters] = double_array_to_struct_array(parameters,'muscle_strength',muscle_strength);
[parameters] = double_array_to_struct_array(parameters,'muscle_pass_stiff_shift',muscle_pass_stiff_shift);
[parameters] = double_array_to_struct_array(parameters,'muscle_pass_stiff_scale',muscle_pass_stiff_scale);
[parameters] = double_array_to_struct_array(parameters,'fiber_damping',fiber_damping);
[parameters] = double_array_to_struct_array(parameters,'muscle_mass',[]);

muscle_info.parameters = parameters;

%% scale muscle-tendon parameters based on user-defined settings
muscle_info = scale_MTparameters(S,muscle_info);

%% calculate muscle-tendon parameter values that depend on others
for i=1:NMuscle
    % shift tendon stiffness curve based on its stiffness
    tendon_stiff_shift_i = getShift(muscle_info.parameters(i).tendon_stiff);
    muscle_info.parameters(i).tendon_stiff_shift = tendon_stiff_shift_i;

    % compute muscle mass
    muscle_mass_i = GetMuscleMass(muscle_info.parameters(i).FMo,...
        muscle_info.parameters(i).lMo,muscle_info.parameters(i).specific_tension);
    muscle_info.parameters(i).muscle_mass = muscle_mass_i;
end


model_info.muscle_info = muscle_info;

% Time constants of excitation-activation dynamics
model_info.muscle_info.tact = S.subject.muscle_activation_time_constant; % activation
model_info.muscle_info.tdeact = S.subject.muscle_deactivation_time_constant; % deactivation

end