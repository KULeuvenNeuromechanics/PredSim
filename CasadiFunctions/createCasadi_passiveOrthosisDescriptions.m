function [f_passiveOrthosis] = createCasadi_passiveOrthosisDescriptions(S,model_info,S_orthosis_settings_i)
% --------------------------------------------------------------------------
% createCasadi_passiveOrthosisDescriptions
%   This function wraps "passiveOrthosisDescriptions" in a CasADi function
%   and returns the handle.
% 
% INPUT:
%   - S -
%   * setting structure S
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
%   - S_orthosis_settings_i -
%   * setting structure for the orthosis. "S_orthosis_settings_i.name" should
%   contain the name of the function you want to call when adding an
%   orthosis to your simulation model. 
%   To add an orthosis, in main.m, set "S.orthosis.settings{i}.name" to the 
%   name of your desired orthosis.
%   Add other fields than "name" to pass settings (e.g. stiffness
%   parameters) to your orthosis description.
%
%
% OUTPUT:
%   - f_passiveOrthosis -
%   * handle of the CasADi function "settings_orthosis.name"
% 
% Original author: Lars D'Hondt
% Original date: 9/September/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.$

n_coord = model_info.ExtFunIO.jointi.nq.all;

%% Define variables
% all inputs
q = SX.sym('q',n_coord);
qdot = SX.sym('qdot',n_coord);
GRFy = SX.sym('GRFy',2);

% all outputs
T = SX('T',n_coord,1);

% indices of used inputs
idx_q = [];
idx_qdot = [];
idx_GRFy = [];

% get function that describes orthosis
fun_passiveOrthosis = passiveOrthosisDescriptions(S_orthosis_settings_i);

% identify function inputs
n_arg_in = nargin(fun_passiveOrthosis);

for j=1:n_arg_in
    


end










