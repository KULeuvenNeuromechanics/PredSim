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
%       Default value = 35, Zajac (1989). This value can be scaled by a
%       factor given in "S.subject.tendon_stiff"
%       - tensions: specific muscle tension. Data from Uchida et al. (2016)
%       is loaded through getSpecificTensions.m
%       - pctsts: fraction of muscle fibers that are slow twitch fibers. Data
%       from Uchida et al. (2016) is loaded through getSlowTwitchRatios.m
%       
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

% This function reads all parameters that describe the muscle- and tendon
% properties. 

muscleNames = model_info.muscle_info.muscle_names;
NMuscle = model_info.muscle_info.NMuscle;

%% read default muscle- and tendon parameters from opensim model file
t0 = tic;
[FMo, lMo, lTs, alphao, vMmax] = getMTparameters(osim_path,muscleNames);
disp(['   reading MT params: ' num2str(toc(t0)) ' s'])

% maximum isometric force (N)
model_info.muscle_info.FMo = FMo;

% optimal fiber length (m)
model_info.muscle_info.lMo = lMo;

% tendon slack length (m)
model_info.muscle_info.lTs = lTs;

% pennation angle at optimal fiber length (rad)
model_info.muscle_info.alphao = alphao;

% maximum contraction velocity (m/s)
model_info.muscle_info.vMmax = vMmax;

%% parameters not from osim file
% default tendon stiffness
model_info.muscle_info.aTendon = 35*ones(NMuscle,1);

% specific tensions of muscle fibers
model_info.muscle_info.tensions = getSpecificTensions(muscleNames);

% ratio of slow twitch muscle fibers
model_info.muscle_info.pctsts = getSlowTwitchRatios(muscleNames);

%% scale muscle-tendon parameters based on user-defined settings
t0 = tic;
model_info = scale_MTparameters(S,model_info);
disp(['   scaling MT params: ' num2str(toc(t0)) ' s'])

%% impose symmetry on the muscle-tendon parameters
if ~isempty(S.subject.muscle_sym) && S.subject.muscle_sym
    t0 = tic;
    model_info = impose_symmetry_MTparameters(S,model_info);
    disp(['   MT params symmetry: ' num2str(toc(t0))])
end

%% calculate muscle-tendon parameter values that depend on others
% shift tendon stiffness curve based on its stiffness
model_info.muscle_info.shift = getShift(model_info.muscle_info.aTendon);

% compute muscle mass
model_info.muscle_info.muscle_mass = GetMuscleMass(model_info.muscle_info.FMo,...
    model_info.muscle_info.lMo, model_info.muscle_info.tensions);






