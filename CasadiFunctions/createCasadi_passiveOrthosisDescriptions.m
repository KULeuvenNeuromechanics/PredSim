function [f_passiveOrthosis] = createCasadi_passiveOrthosisDescriptions(S,model_info,S_orthosis_settings_i)
% --------------------------------------------------------------------------
% createCasadi_passiveOrthosisDescriptions
%   This function wraps a selected function from "passiveOrthosisDescriptions" 
%   in a CasADi function and returns the handle.
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
%   Edit: Added support for GRF input to orthosis
% Last edit by: Lars D'Hondt
% Last edit date: 19/October/2022
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
GRF = SX.sym('GRF',6);

% all outputs
T = SX(n_coord,1);

% selected inputs
orth_in = cell(1,n_arg_in);
orth_in{1} = S_orthosis_settings_i;
T_out = cell(1,n_arg_out);

%% Identify inputs
% loop over inputs
for i=1:length(names_in)
    name_i = names_in{i};
    if strcmp(name_i(1:5),'qdot_')
        % index of coordinate that corresponds to asked input
        idx = find(strcmp(coord_names,{name_i(6:end)}));
        if length(idx)==1
            % add velocity of coordinate to orthosis inputs
            orth_in{i+1} = qdot(idx);
        else
            error(['Invalid description of orthosis "' oname,...
                '". Coordinate "' name_i(6:end) '" not uniquely defined in OpenSim model.'])
        end

    elseif strcmp(name_i(1:2),'q_')
        % index of coordinate that corresponds to asked input
        idx = find(strcmp(coord_names,{name_i(3:end)}));
        if length(idx)==1
            % add velocity of coordinate to orthosis inputs
            orth_in{i+1} = q(idx);
        else
            error(['Invalid description of orthosis "' oname,...
                '". Coordinate "' name_i(3:end) '" not uniquely defined in OpenSim model.'])
        end

    elseif strcmp(name_i(1:3),'GRF')
        % index of ground reaction force vector component that corresponds to asked input
        idx_c = strfind('xyz',name_i(5));

        if isempty(idx)
            error(['Invalid description of orthosis "' oname,...
                '". Please use "GRF_*_r" or "GRF_*_l" (where * is "x", "y", or "z") ',...
                'to refer to a component of the total ground reaction force on a foot.'])
        end
        
        % GRF from left or right foot
        idx_lr = strfind('lr',name_i(7));
        if isempty(idx_lr)
            error(['Invalid description of orthosis "' oname,...
                '". Please use e.g. "GRF_y_r" or "GRF_y_l" to refer to a the vertical component of ',...
                'the total ground reaction force on the right or left foot respectively.'])
        end

        % GRF order is left, right
        idx = idx_c + 3*(idx_lr-1);
        % add GRF component to orthosisinputs
        orth_in{i+1} = GRF(idx);

    else
        error(['Invalid description of orthosis "' oname,...
                '". Input names should start with "q_" or "qdot_", followed by coordinate name.'])
    end
end % end loop over inputs

%% Identify outputs
idx_T = [];
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
[T_out{:}] = fun_passiveOrthosis(orth_in{:});

% assign outputs
for i=1:n_arg_out-2
    T(idx_T(i)) = T_out{i+2};
end


%% Create CasADi function

f_passiveOrthosis = Function('f_passiveOrthosis',{q,qdot,GRF},{T},{'q','qdot','GRF'},{'T'});










