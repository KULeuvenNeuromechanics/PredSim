function [model_info] = get_passive_moment_info(S,model_info)
% --------------------------------------------------------------------------
% get_passive_moment_info
%   This function organises the different coefficients that describe
%   passive torques and moments. They are added to model_info, to be used
%   when creating the CasADi functions
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
% Original author: Lars D'Hondt
% Original date: 15/April/2022
%
% Last edit by:
% Last edit date:
% --------------------------------------------------------------------------

% number of coordinates
n_coord = model_info.ExtFunIO.jointi.nq.all;
% coordinate names
coord_names = model_info.ExtFunIO.coord_names.all;

%% Collect all coefficients
% damping coefficients
if ~isempty(S.subject.set_damping_coefficient_selected_dofs)
    try
        damping_coeff = unpack_name_value_combinations(S.subject.set_damping_coefficient_selected_dofs,...
            coord_names,1);
    catch errmsg
        error(['Unable to set damping coefficient because: ', errmsg]);
    end
end

% stiffness coefficients
if ~isempty(S.subject.set_stiffness_coefficient_selected_dofs)
    try
        stiffness_coeff = unpack_name_value_combinations(S.subject.set_stiffness_coefficient_selected_dofs,...
            coord_names,1);
     catch errmsg
        error(['Unable to set stiffness coefficient because: ', errmsg.message]);
    end
end

% limit torque coefficients
[K_pass_default,theta_pass_default] = get_default_coord_limit_torque_coefficients(coord_names);
if ~isempty(S.subject.set_limit_torque_coefficients_selected_dofs)
    try
        [K_pass,theta_pass] = unpack_name_value_combinations(S.subject.set_limit_torque_coefficients_selected_dofs,...
            coord_names,[4,2]);
        K_pass = K_pass';
        theta_pass = theta_pass';
     catch errmsg
        error(['Unable to set limit torque coefficients because: ', errmsg.message]);
    end
end


%% Organise in struct
limitTorque = [];

for i=1:n_coord
    if max(model_info.ExtFunIO.jointi.floating_base == i)>0
        passive_moment_info(i).coord_name = coord_names{i};
        passive_moment_info(i).K_pass = [];
        passive_moment_info(i).theta_pass = [];
        passive_moment_info(i).damping_coeff = 0;
        passive_moment_info(i).stiffness_coeff = 0;

    else
        passive_moment_info(i).coord_name = coord_names{i};
        if exist("K_pass","var") && ~isnan(K_pass(i,1))
            passive_moment_info(i).K_pass = K_pass(i,:);
            passive_moment_info(i).theta_pass = theta_pass(i,:);
        elseif exist("K_pass","var") && K_pass(i,1)==0
            passive_moment_info(i).K_pass = [];
            passive_moment_info(i).theta_pass = [];
        elseif ~isnan(K_pass_default(i,1))
            passive_moment_info(i).K_pass = K_pass_default(i,:);
            passive_moment_info(i).theta_pass = theta_pass_default(i,:);
        else
            passive_moment_info(i).K_pass = [];
            passive_moment_info(i).theta_pass = [];
        end

        if exist("damping_coeff","var") && ~isnan(damping_coeff(i))
            passive_moment_info(i).damping_coeff = damping_coeff(i);
        else
            passive_moment_info(i).damping_coeff = S.subject.damping_coefficient_all_dofs;
        end

        if exist("stiffness_coeff","var") && ~isnan(stiffness_coeff(i))
            passive_moment_info(i).stiffness_coeff = stiffness_coeff(i);
        else
            passive_moment_info(i).stiffness_coeff = S.subject.stiffness_coefficient_all_dofs;
        end

    end

    if ~isempty(passive_moment_info(i).K_pass) && ~isempty(passive_moment_info(i).theta_pass)
        limitTorque(end+1) = i;
    end
end

model_info.passive_moment_info.parameters = passive_moment_info;


% Coordinates with limit torque
model_info.ExtFunIO.jointi.limitTorque = limitTorque;
model_info.ExtFunIO.jointi.nq.limTorq = length(limitTorque);


end