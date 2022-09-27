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

import casadi.*

n_coord = model_info.ExtFunIO.jointi.nq.all;
coord_names = model_info.ExtFunIO.coord_names.all;

%% Check provided function
% orthosis name
oname = S_orthosis_settings_i.name;

% get function that describes orthosis
fun_passiveOrthosis = passiveOrthosisDescriptions(S_orthosis_settings_i);

% size of input and output
n_arg_in = nargin(fun_passiveOrthosis);
n_arg_out = nargout(fun_passiveOrthosis);

% dummy input
dummy_input = cell(1,n_arg_in);
dummy_input{1} = S_orthosis_settings_i;
dummy_input(2:end) = {0};

% dummy output
dummy_output = cell(1,n_arg_out);

% evaluate fun with dummy
[dummy_output{:}] = fun_passiveOrthosis(dummy_input{:});

% restructure outputs
names_in = dummy_output{1};
names_out = dummy_output{2};
torques = dummy_output(3:end);

% check input names
if length(names_in) ~= n_arg_in-1
    error(['Invalid description of orthosis "' oname,...
        '". Number of input names does not match number of inputs.'])
end

% check output names
if length(names_out) ~= length(torques)
    error(['Invalid description of orthosis "' oname,...
        '". Number of output names does not match number of outputs.'])
end

%% Define variables
% all inputs
q = SX.sym('q',n_coord);
qdot = SX.sym('qdot',n_coord);

% all outputs
T = SX(n_coord,1);

% selected inputs
q_qdot_in = cell(1,n_arg_in);
q_qdot_in{1} = S_orthosis_settings_i;
T_out = cell(1,n_arg_out);

%% Identify inputs and outputs
% indices of used inputs and outputs
idx_q = [];
idx_qdot = [];
idx_T = [];

% identify inputs
for i=1:length(names_in)
    name_i = names_in{i};
    if strcmp(name_i(1:5),'qdot_')
        idx = find(strcmp(coord_names,{name_i(6:end)}));
        if length(idx)==1
            idx_qdot(end+1) = idx;
            q_qdot_in{i+1} = qdot(idx);
        else
            error(['Invalid description of orthosis "' oname,...
                '". Coordinate "' name_i(6:end) '" not uniquely defined in OpenSim model.'])
        end
    elseif strcmp(name_i(1:2),'q_')
        idx = find(strcmp(coord_names,{name_i(3:end)}));
        if length(idx)==1
            idx_q(end+1) = idx;
            q_qdot_in{i+1} = q(idx);
        else
            error(['Invalid description of orthosis "' oname,...
                '". Coordinate "' name_i(3:end) '" not uniquely defined in OpenSim model.'])
        end
    else
        error(['Invalid description of orthosis "' oname,...
                '". Input names should start with "q_" or "qdot_", followed by coordinate name.'])
    end
end

% identify outputs
for i=1:length(names_out)
    name_i = names_out{i};
    idx = find(strcmp(coord_names,{name_i}));
    if length(idx)==1
        idx_T(end+1) = idx;
    else
        error(['Invalid description of orthosis "' oname,...
            '". Coordinate "' name_i '" not uniquely defined in OpenSim model.'])
    end
end


%% Evaluate function
[T_out{:}] = fun_passiveOrthosis(q_qdot_in{:});

% assign outputs
for i=1:n_arg_out-2
    T(idx_T(i)) = T_out{i+2};
end


%% Create CasADi function

f_passiveOrthosis = Function('f_passiveOrthosis',{q,qdot},{T},{'q','qdot'},{'T'});










