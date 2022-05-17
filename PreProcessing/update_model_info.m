function [model_info] = update_model_info(S,osim_path,model_info)
% --------------------------------------------------------------------------
% update_model_info
%   Create additional fields with indices combination in model_info, based 
%   on existing fields. All sets of indices needed outside PreProcessing 
%   should be added here to guarantee consistency.
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% OUTPUT:
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% Original author: Dhruv Gupta
% Original date: 17/March/2022
%
% Last update:
%   moved most content to separate functions that are called earlier
% Last edit by: Lars D'Hondt
% Last edit date: 12/April/2022
% --------------------------------------------------------------------------


%% Add crossed coordinates to muscle_info
% This is purely for user convenience, functions use muscle_spanning_joint_info

NMuscle = model_info.muscle_info.NMuscle;
for i=1:NMuscle
    idx_coords_i = find(model_info.muscle_info.muscle_spanning_joint_info(i,:)==1);
    names_coords_i = model_info.ExtFunIO.coord_names.all(idx_coords_i);
    tmpst = [num2str(length(names_coords_i)) ': '];
    for j=1:length(names_coords_i)-1
        tmpst = [tmpst names_coords_i{j} ', '];
    end
    tmpst = [tmpst names_coords_i{end}];
    model_info.muscle_info.parameters(i).coords_crossed_by_muscle = tmpst;
end


%% Coordinate index subsets
Ncoords = length(model_info.ExtFunIO.coord_names.all);

model_info.ExtFunIO.jointi.armsi = sort([model_info.ExtFunIO.jointi.arm_l model_info.ExtFunIO.jointi.arm_r]);
model_info.ExtFunIO.jointi.noarmsi = setdiff(1:Ncoords,model_info.ExtFunIO.jointi.armsi);
model_info.ExtFunIO.jointi.base_forward = model_info.ExtFunIO.symQs.base_forward;
model_info.ExtFunIO.jointi.base_lateral = model_info.ExtFunIO.symQs.base_lateral;
% model_info.ExtFunIO.jointi.base_forward = setdiff(model_info.ExtFunIO.symQs.QdotsInvA,model_info.ExtFunIO.symQs.QsInvA);

model_info = addCoordNames(model_info,'muscleActuated');

%% Number of degrees of freedom for later use
model_info.ExtFunIO.jointi.nq.all          = Ncoords; % all
model_info.ExtFunIO.jointi.nq.arms         = length(model_info.ExtFunIO.jointi.armsi);
model_info.ExtFunIO.jointi.nq.noArms       = length(model_info.ExtFunIO.jointi.noarmsi);
model_info.ExtFunIO.jointi.nq.musAct       = length(model_info.ExtFunIO.jointi.muscleActuated);
model_info.ExtFunIO.jointi.nq.torqAct      = length(model_info.ExtFunIO.jointi.torqueActuated);
model_info.ExtFunIO.jointi.nq.rot          = length(model_info.ExtFunIO.jointi.rotations);
model_info.ExtFunIO.jointi.nq.trnsl        = length(model_info.ExtFunIO.jointi.translations);

%% Model symmetry (for half gait cycle simulations)
ActOpp = find(ismember(model_info.ExtFunIO.symQs.QsOpp(:),model_info.ExtFunIO.jointi.torqueActuated));
ActInvA = setdiff(1:model_info.ExtFunIO.jointi.nq.torqAct,ActOpp);

ActInvB = zeros(size(ActInvA));

for i=1:length(ActInvB)
    idx_i = find(model_info.ExtFunIO.symQs.QsInvA==model_info.ExtFunIO.jointi.torqueActuated(ActInvA(i)));
    ActInvB(i) = find(model_info.ExtFunIO.jointi.torqueActuated==model_info.ExtFunIO.symQs.QsInvB(idx_i));
end
model_info.ExtFunIO.symQs.ActInvA = ActInvA;
model_info.ExtFunIO.symQs.ActInvB = ActInvB;
model_info.ExtFunIO.symQs.ActOpp = ActOpp;


%% Pelvis height used for initial guess
model_info = get_IG_pelvis_y(S,osim_path,model_info);

% add osim_path so it will be included in saved results
model_info.osim_path = osim_path;










