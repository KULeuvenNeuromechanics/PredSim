function [muscle_info] = scale_MTparameters(S,muscle_info)
% --------------------------------------------------------------------------
% scale_MTparameters
%   This functions scales the muscle-tendon parameter values. See
%   "read_and_scale_MTparameters" for information on these parameters.
%   Scaling depends on the user settings.
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - muscle_info -
%   * structure with all the parameters descibing muscle-tendon dynamics
% 
% OUTPUT:
%   - muscle_info -
%   * structure with all the parameters descibing muscle-tendon dynamics
% 
% Original author: Lars D'Hondt
% Original date: 17/March/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

muscleNames = muscle_info.muscle_names;
NMuscle = muscle_info.NMuscle;

%% Unpack scaling settings of muscle-tendon parameters from .osim file
% The scale factors for FMo, lMo, lTs, alphao, and vMmax are input as 
% "S.subject.scale_MT_params = {{'soleus_l'},'FMo',0.9,{'soleus_l'},'alphao',1.1};"
% We reshape this to "scale_MTparams.FMo = {{'soleus_l'},0.9}" and 
% "scale_MTparams.alphao = {{'soleus_l'},1.1}".

if ~isempty(S.subject.scale_MT_params)
    try
        % Test length of input
        if mod(length(S.subject.scale_MT_params),3)~=0
            error('Expected input pattern is: property, muscle names, scale factor')
        end
        % Allocate scale_MTparams
        scale_MTparams.FMo = {};
        scale_MTparams.lMo = {};
        scale_MTparams.lTs = {};
        scale_MTparams.alphao = {};
        scale_MTparams.vMmax = {};
        MTproperties = fieldnames(scale_MTparams);
        % Loop over sets of 3 inputs
        for i=1:length(S.subject.scale_MT_params)/3
            i_mus = (i-1)*3+1; % index of muscle names
            i_prp = i_mus+1; % index of property
            i_scf = i_mus+2; % index of scale factor
            % Test property is valid
            idx_property = find(strcmp(MTproperties(:),S.subject.scale_MT_params(i_prp)));
            if ~isempty(idx_property)
                % Place muscle names and scale factor in field of
                % scale_MTparams that corresponds to property to be scaled
                scale_MTparams.(MTproperties{idx_property}){end+1} = S.subject.scale_MT_params{i_mus};
                scale_MTparams.(MTproperties{idx_property}){end+1} = S.subject.scale_MT_params{i_scf};
            else
                % Throw error if property is invalid
                error([S.subject.scale_MT_params{i_prp} ' is not an accepted muscle-tendon property. ',...
                    'Possible entries are: FMo, lMo, lTs, alphao, vMmax.'])
            end
        end
    catch errmsg
        error(['Unable to scale muscle-tendon parameters because: ', errmsg.message]);
    end
end

%% Add scaling settings of non-osim muscle-tendon parameters
% muscle strength, passive muscle stiffness, and tendon stiffness are
% parameters related to the implementation of the muscle model, and are not
% uniquely related to muscle-tendon parameters found in an .osim file.
% Each of their scale factors is defined as a separate setting. Here, the
% scale factors are added to the struct defined in the previous section.
scale_MTparams.muscle_strength = S.subject.muscle_strength;
scale_MTparams.muscle_pass_stiff_scale = S.subject.muscle_pass_stiff_scale;
scale_MTparams.muscle_pass_stiff_shift = S.subject.muscle_pass_stiff_shift;
scale_MTparams.tendon_stiff = S.subject.tendon_stiff_scale;


%% Scale all muscle-tendon parameters

% Loop over all parameters
MTproperties = fieldnames(scale_MTparams);
for i=1:length(MTproperties)
    % Test the field is empty i.e. this parameter does not need to be
    % scaled for any muscle
    if ~isempty(scale_MTparams.(MTproperties{i}))
        try
            % Get array with scale factors for every muscle
            scale_factors_MTparam_i = get_scale_factors(scale_MTparams.(MTproperties{i}), NMuscle, muscleNames);
            % Multiply parameter with scale factor
            for j=1:NMuscle
                MTparam_ij = muscle_info.parameters(j).(MTproperties{i})*scale_factors_MTparam_i(j);
                muscle_info.parameters(j).(MTproperties{i}) = MTparam_ij;
            end
        catch errmsg
            error(['Unable to scale muscle-tendon parameter "' MTproperties{i} '" because: ', errmsg.message]);
        end
    end
end


%% helper functions
% get array of scale factors
function [scale_factors] = get_scale_factors(S_field, N_muscle, muscle_names)
    % Scale factors are given as a cell array with alternating musclenames 
    % and scale factor. Call helper from /VariousFunctions to rewrite this
    % as an array of number-of-muscles doubles, with the scale factors
    % placed in entries corresponding to the defined order of the muscles.
    scale_factors = unpack_name_value_combinations(S_field, muscle_names, 1);
    
    % This helper function has nan as default output values. Replace this
    % by scale factor 1.
    for ii=1:N_muscle
        if isnan(scale_factors(ii))
            scale_factors(ii) = 1;
        else
            scale_factors(ii) = to_scale_factor(scale_factors(ii));
        end
    end

end
% transform given scale factor to valid scale factor
function[sf] = to_scale_factor(some_number)
    % Scale factor is always positive
    if some_number >=0
        sf = some_number;
    else
        error(['Scale factor ' num2str(some_number) ' should be positive.'])
    end

    % Might want to add check that scale factor is not excessively large...
end
end