function [model_info] = get_model_info(S,osim_path)
% --------------------------------------------------------------------------
% get_model_info
%   Create model_info, a struct that contains all information about the
%   model. The indices of joints and coordinates are loaded from *_IO.mat, 
%   which was created by osim2dll. The muscle names are read directly 
%   from the provided opensim model. 
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
%
% OUTPUT:
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% Original author: Lars D'Hondt
% Original date: 11/April/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


% extract name of opensim model file
[~,osim_file_name,~] = fileparts(osim_path);

%% Subject mass
if isempty(S.subject.mass)
    model_info.mass = getModelMass(osim_path);
else
    model_info.mass = S.subject.mass;
end

%% External function
% load IO from external function
% We load IO into a struct, because having a variable named "IO" conflicts
% with OpenSim api.
IEAIAIO = load(fullfile(S.misc.subject_path,['F_' osim_file_name '_IO.mat']),'IO');
ExtFunIO = IEAIAIO.IO;

% create model_info with IO inside
model_info.ExtFunIO = ExtFunIO;

% coordinate names
model_info.ExtFunIO.coord_names.all = fieldnames(model_info.ExtFunIO.coordi);

% number of coordinates
model_info.ExtFunIO.jointi.nq.all = length(model_info.ExtFunIO.coord_names.all);

% all inputs
fields = ["Qs", "Qdots", "Qdotdots"];
for i=fields
    for j=1:length(model_info.ExtFunIO.coord_names.all)
        model_info.ExtFunIO.input.(i).all(j) =...
            model_info.ExtFunIO.input.(i).(model_info.ExtFunIO.coord_names.all{j});
    end
end

%% OpenSim model file
% read muscle names from .osim file
import org.opensim.modeling.*;
model = Model(osim_path);

muscle_names = cell(1,model.getMuscles().getSize());
for i=1:model.getMuscles().getSize()
    muscle_names{i} = char(model.getMuscles().get(i-1).getName());
end

model_info.muscle_info.muscle_names = muscle_names;

model_info.muscle_info.NMuscle = length(muscle_names);

%% OpenSim API
% indices of coordinates in the OpenSim API state vector
model_info = getCoordinateIndexForStateVectorOpenSimAPI(S,osim_path,model_info);

%% Symmetry
[symQs, model_info.ExtFunIO.jointi] = identify_kinematic_chains(S,osim_path,model_info);

orderMus = 1:length(model_info.muscle_info.muscle_names);
orderMusInv = zeros(1,length(model_info.muscle_info.muscle_names));
for i=1:length(model_info.muscle_info.muscle_names)

    orderMusInv(i) = find(strcmp(model_info.muscle_info.muscle_names,...
        mirrorName(model_info.muscle_info.muscle_names{i})));
end
symQs.MusInvA = orderMus;
symQs.MusInvB = orderMusInv;

model_info.ExtFunIO.symQs = symQs;





end % end of function
