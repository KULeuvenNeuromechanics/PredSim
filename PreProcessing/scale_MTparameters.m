function [model_info] = scale_MTparameters(S,model_info)
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
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% OUTPUT:
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% Original author: Lars D'Hondt
% Original date: 17/March/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

muscleNames = model_info.muscle_info.muscle_names;
NMuscle = model_info.muscle_info.NMuscle;

%% scaling muscle strength
if ~isempty(S.subject.muscle_strength)
    try
        % get scale factor array
        scale_factors_muscle_strength = get_scale_factors(S.subject.muscle_strength, NMuscle, muscleNames);
        % scale
        model_info.muscle_info.muscle_strength = model_info.muscle_info.muscle_strength.*scale_factors_muscle_strength;
    catch errmsg
        error(['Unable to scale muscle strength because: ', errmsg]);
    end
end

%% scaling muscle stiffness
if ~isempty(S.subject.muscle_stiff)
    try
        % get scale factor array
        scale_factors_muscle_stiff = get_scale_factors(S.subject.muscle_stiff, NMuscle, muscleNames);
        % scale
        model_info.muscle_info.muscle_stiffness = model_info.muscle_info.muscle_stiffness.*scale_factors_muscle_stiff;
    catch errmsg
        error(['Unable to scale muscle stiffness because: ', errmsg]);
    end
end

%% scaling tendon stiffness
if ~isempty(S.subject.tendon_stiff)
    try
        % get scale factor array
        scale_factors_tendon_stiff = get_scale_factors(S.subject.tendon_stiff, NMuscle, muscleNames);
        % scale
        model_info.muscle_info.aTendon = model_info.muscle_info.aTendon.*scale_factors_tendon_stiff;
    catch errmsg
        error(['Unable to scale tendon stiffness because: ', errmsg]);
    end
end

%% scaling muscle-tendon parameters from opensim model file
if ~isempty(S.subject.MT_params)
    try
        if mod(length(S.subject.MT_params),3)~=0
            error('Expected input pattern is: property, muscle names, scale factor')
        end
        scale_MTparams.FMo = {};
        scale_MTparams.lMo = {};
        scale_MTparams.lTs = {};
        scale_MTparams.alphao = {};
        scale_MTparams.vMmax = {};
        MTproperties = fieldnames(scale_MTparams);
        for i=1:mod(length(S.subject.MT_params),3)
            idx_property = find(strcmp(MTproperties{:},S.subject.MT_params{i}));
            if ~isempty(idx_property)
                scale_MTparams.(MTproperties{idx_property}){end+1} = S.subject.MT_params{i+1};
                scale_MTparams.(MTproperties{idx_property}){end+1} = S.subject.MT_params{i+2};
            else
                error([S.subject.MT_params{i} 'is not an accepted muscle-tendon property. ',...
                    'Possible entries are: FMo, lMo, lTs, alphao, vMmax.'])
            end
        end
        for i=1:5
            % get scale factor array
            scale_factors_MTparam_i = get_scale_factors(scale_MTparams.(MTproperties{i}), NMuscle, muscleNames);
            % scale
            model_info.muscle_info.(MTproperties{i}) = model_info.muscle_info.(MTproperties{i}).*scale_factors_MTparam_i;
        end
    catch errmsg
        error(['Unable to scale muscle-tendon parameters because: ', errmsg]);
    end
end

%% helper functions
% get array of scale factors
function [scale_factors] = get_scale_factors(S_field, N_muscle, muscle_names)

    scale_factors = unpack_name_value_combinations(S_field, muscle_names, 1);
    
    for ii=1:N_muscle
        if isnan(scale_factors(ii))
            scale_factors(ii) = 1;
        else
            scale_factors(ii) = to_scale_factor(scale_factors(ii));
        end
    end

% N_fields = length(S_field);
% if mod(N_fields,2)~=0
%     error('Array of name-value pairs of muscle name and scale factor needs to have an even number of inputs.')
% end
% scale_factors = ones(1,N_muscle);
% 
% for ii=1:N_fields/2
%     muscle_names_i = S_field{2*ii-1};
%     if ~iscell(muscle_names_i)
%         muscle_names_i = {muscle_names_i};
%     end
%     for jj=1:length(muscle_names_i)
%         idx_muscle_ij = find(contains(muscle_names{:},muscle_names_i{jj}));
%         if isempty(idx_muscle_ij)
%             error(['Muscle name not found: ' muscle_names_i{jj}])
%         end
%         if scale_factors(idx_muscle_ij)==1
%             scale_factors(idx_muscle_ij) = to_scale_factor(S_field{2*ii});
%         else
%             error(['Cannot scale ' muscle_names{idx_muscle_ij} ' twice.'])
%         end
%     end
% end

end
% transform given scale factor to valid scale factor
function[sf] = to_scale_factor(some_number)
    if some_number >=0 && some_number <= 1
        sf = some_number;
    elseif some_number >=0 && some_number <= 100
        sf = some_number/100;
        warning(['Assuming scale factor input value as ' num2str(some_number),...
            ' to be a percentage. Converting to scale factor = ' num2str(sf) '.'])
    else
        error(['Scale factor ' num2str(some_number) ' is outside the expected range of 0 to 1.'])
    end
end
end