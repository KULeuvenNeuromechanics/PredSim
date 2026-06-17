function [model_info] = getCoordinateIndexForStateVectorOpenSimAPI(S,osim_path,model_info)
% --------------------------------------------------------------------------
% getCoordinateIndexForStateVectorOpenSimAPI
%   This functions finds the index of each coordinate position in the state
%   vector used by the OpenSim API. 
%   The indices are added to model_info: model_info.ExtFunIO.coordi_OpenSimAPIstate
%   This does not contain odd numbers, because those are used for the
%   velocity states.
%   Do note that the OpenSim API starts indexing at 0. 
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
% Original date: 06/April/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

% number of coordinates
n_coord = model_info.ExtFunIO.jointi.nq.all;

%% Initialise model
import org.opensim.modeling.*;
model = Model(osim_path);
state = model.initSystem;
% Get state vector
state_vars = model.getStateVariableValues(state);
for i=1:n_coord
    state_vars.set((i-1)*2,(i-1)*2);
end
model.setStateVariableValues(state,state_vars);
model.realizePosition(state);

%%
coordset = model.getCoordinateSet();

for i=1:n_coord
    coordname_i = model_info.ExtFunIO.coord_names.all{i};
    coord_i = coordset.get(coordname_i);
    coord_state_idx = coord_i.getStateVariableValues(state).getAsMat;
    
    coordi_OpenSimAPI.(coordname_i) = coord_state_idx(1);
end

model_info.ExtFunIO.coordi_OpenSimAPIstate = coordi_OpenSimAPI;

end