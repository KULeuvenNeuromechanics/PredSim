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


%% Read from osim model
import org.opensim.modeling.*;

model = Model(osim_path);
model.initSystem;


PCSA_force = nan(1,model_info.ligament_info.NLigament);
slack_length = PCSA_force;

% loop over all ligaments
for i=1:model_info.ligament_info.NLigament
    lig = model.getForceSet().get(model_info.ligament_info.ligament_names{i});
    lig = Ligament.safeDownCast(lig);
    PCSA_force(1,i) = lig.get_pcsa_force();
    slack_length(1,i) = lig.getRestingLength();

end

%
[parameters] = double_array_to_struct_array([],'ligament_name',model_info.ligament_info.ligament_names);
[parameters] = double_array_to_struct_array(parameters,'cross_section_area',PCSA_force);
[parameters] = double_array_to_struct_array(parameters,'slack_length',slack_length);


%% Set force-length functions
% selected ligaments from settings
ligament_stiffness = unpack_name_value_combinations(S.subject.set_stiffness_selected_ligaments,...
    model_info.ligament_info.ligament_names, 1);

% use default for all other
ligament_stiffness(ismissing(ligament_stiffness)) = S.subject.stiffness_all_ligaments;

[parameters] = double_array_to_struct_array(parameters,'stiffness',ligament_stiffness);

%% Add to model_info
model_info.ligament_info.parameters = parameters;

end