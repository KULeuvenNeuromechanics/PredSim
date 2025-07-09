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
% update:
%   read out ligament names
%
% Last edit by: Lars D'Hondt
% Last edit date: 5/April/2023
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

ligament_names = {};
for i=1:model.getForceSet().getSize()
    force_i = model.getForceSet().get(i-1).getConcreteClassName();
    if strcmp(force_i,'Ligament')
        ligament_names{1,end+1} = char(model.getForceSet().get(i-1).getName());
    end
end
model_info.ligament_info.ligament_names = ligament_names;
model_info.ligament_info.NLigament = length(ligament_names);

%% OpenSim API
% indices of coordinates in the OpenSim API state vector
model_info = getCoordinateIndexForStateVectorOpenSimAPI(S,osim_path,model_info);

%% Symmetry
[symQs, model_info.ExtFunIO.jointi] = identify_kinematic_chains(S,osim_path,model_info);

orderMus = 1:length(model_info.muscle_info.muscle_names);
orderMusInv = orderMus;
for i=1:length(model_info.muscle_info.muscle_names)

    idx_mus_inv_i = find(strcmp(model_info.muscle_info.muscle_names,...
        mirrorName(model_info.muscle_info.muscle_names{i})));
    if ~isempty(idx_mus_inv_i)
        orderMusInv(i) = idx_mus_inv_i;
    else
        if strcmpi(S.misc.gaitmotion_type,'HalfGaitCycle')
            error("Model asymmetry detected while S.misc.gaitmotion_type = " + ...
                "'HalfGaitCycle'")
        end
    end
end
symQs.MusInvA = orderMus;
symQs.MusInvB = orderMusInv;

model_info.ExtFunIO.symQs = symQs;

%% indices for left and right side muscles
idx_mus_l = [];
idx_mus_r = [];

for i=1:length(model_info.muscle_info.muscle_names)

    mus_i = model_info.muscle_info.muscle_names{i};
    [mus_i_l, mus_i_r] = mirrorName(mus_i);
    if strcmp(mus_i, mus_i_r)
        idx_mus_r(end+1) = i;

        % For half GC, set idx of left muscles to be symmetric to idx of
        % right muscles.
        if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
            idx_mus_l(end+1) = find(strcmp(model_info.muscle_info.muscle_names, mus_i_l));

        end

        % For full GC, we do not need the idx to be symmetric. Additionally,
        % some muscle might only exist left or right.
    elseif ~strcmp(S.misc.gaitmotion_type,'HalfGaitCycle') && strcmp(mus_i, mus_i_l)
        idx_mus_l(end+1) = i;

    end
end

model_info.muscle_info.idx_left = idx_mus_l;
model_info.muscle_info.idx_right = idx_mus_r;


%% add osim_path so it will be included in saved results
model_info.osim_path = osim_path;



end % end of function
