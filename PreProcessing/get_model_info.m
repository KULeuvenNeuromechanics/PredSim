function model_info = get_model_info(S,osim_path)
%
% to do
%
% Author: Lars D'Hondt
%
% Date: 18/January/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,osim_file_name,~] = fileparts(osim_path);

% load IO from external function
load(fullfile(S.misc.subject_folder,['F_' osim_file_name '_IO.mat']),'IO');

model_info.ExtFunIO = IO;

% read muscle names from .osim file
import org.opensim.modeling.*;
model = Model(osim_path);

muscle_names = cell(model.getMuscles().getSize());
for i=1:model.getMuscles().getSize()
    
    muscle_names{i} = char(model.getMuscles().get(i-1).getName());
end

model_info.muscle_info.muscle_names = muscle_names;
