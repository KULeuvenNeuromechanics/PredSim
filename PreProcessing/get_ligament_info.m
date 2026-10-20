function [model_info] = get_ligament_info(S,osim_path,model_info)
% --------------------------------------------------------------------------
% get_ligament_info
%   Read PCSA force and slack length of the ligaments from the osim file. 
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
% Original date: 4/April/2023
%
% Last edit by: Bryce A Killen
% Last edit date: 
% --------------------------------------------------------------------------

% changing this to allow opensim-jam blankevoort ligament type
% https://simtk.org/api_docs/opensim/api_docs/classOpenSim_1_1Blankevoort1991Ligament.html

%% Read from osim model
import org.opensim.modeling.*;

model = Model(osim_path);
model.initSystem;

% set empties 
linear_stiffness = nan(1,model_info.ligament_info.NLigament);
slack_length = linear_stiffness;

% loop over all ligaments
for i=1:model_info.ligament_info.NLigament
    lig = model.getForceSet().get(model_info.ligament_info.ligament_names{i});
    lig = Blankevoort1991Ligament.safeDownCast(lig);
    linear_stiffness(1,i) = lig.get_linear_stiffness();
    slack_length(1,i) = lig.get_slack_length();
    % kind of messed up name - but keeping beucase it is used int he casadi
    % functions/codes
    ligament_stiffness(1,i) = {S.subject.stiffness_all_ligaments};
end

%
[parameters] = double_array_to_struct_array([],'ligament_name',model_info.ligament_info.ligament_names);
[parameters] = double_array_to_struct_array(parameters,'linear_stiffness',linear_stiffness);
[parameters] = double_array_to_struct_array(parameters,'slack_length',slack_length);


%% Set force-length functions
% selected ligaments from settings
%ligament_stiffness = unpack_name_value_combinations(S.subject.set_stiffness_selected_ligaments,...
%    model_info.ligament_info.ligament_names, 1);

% use default for all other
%ligament_stiffness(ismissing(ligament_stiffness)) = S.subject.stiffness_all_ligaments;

[parameters] = double_array_to_struct_array(parameters,'stiffness',ligament_stiffness);

%% Add to model_info
model_info.ligament_info.parameters = parameters;

end