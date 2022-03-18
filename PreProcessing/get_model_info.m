function model_info = get_model_info(S,osim_path)
%
% Create model_info, a struct that contains all information about the
% model.T The muscle properties are read directly from the provided opensim 
% model. The indices of joints and coordinates are loaded from *_IO.mat, 
% which was created by osim2dll.
%
%
% Author: Lars D'Hondt
%
% Date: 18/January/2022
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract name of opensim model file
[~,osim_file_name,~] = fileparts(osim_path);

% load IO from external function
load(fullfile(S.misc.subject_path,['F_' osim_file_name '_IO.mat']),'IO');

% convert indices int32 to doubles
IOfields = fields(IO);
for i=1:length(IOfields)
    IO.(IOfields{i}) = convert2double(IO.(IOfields{i}));
end

% create model_info with IO inside
model_info.ExtFunIO = IO;


% read muscle names from .osim file
import org.opensim.modeling.*;
model = Model(osim_path);

muscle_names = cell(model.getMuscles().getSize());
for i=1:model.getMuscles().getSize()
    muscle_names{i} = char(model.getMuscles().get(i-1).getName());
end

model_info.muscle_info.muscle_names = muscle_names;

model_info.muscle_info.NMuscle = length(muscle_names);


