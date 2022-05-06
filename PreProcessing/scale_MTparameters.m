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

%% unpack scaling settings of muscle-tendon parameters from .osim file
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
    catch errmsg
        error(['Unable to scale muscle-tendon parameters because: ', errmsg]);
    end
end

%% add other scaling settings
scale_MTparams.muscle_strength = S.subject.muscle_strength;
scale_MTparams.muscle_pass_stiff_scale = S.subject.muscle_pass_stiff_scale;
scale_MTparams.muscle_pass_stiff_shift = S.subject.muscle_pass_stiff_shift;
scale_MTparams.tendon_stiff = S.subject.tendon_stiff;


%%
MTproperties = fieldnames(scale_MTparams);
for i=length(MTproperties)
    if ~isempty(scale_MTparams.(MTproperties{i}))
        try
            % get scale factor array
            scale_factors_MTparam_i = get_scale_factors(scale_MTparams.(MTproperties{i}), NMuscle, muscleNames);
            % scale
            for j=1:NMuscle
                MTparam_ij = muscle_info.parameters(j).(MTproperties{i})*scale_factors_MTparam_i(j);
                muscle_info.parameters(j).(MTproperties{i}) = MTparam_ij;
            end
        catch errmsg
            error(['Unable to scale muscle-tendon parameter "' MTproperties{i} '" because: ', errmsg]);
        end
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

end
% transform given scale factor to valid scale factor
function[sf] = to_scale_factor(some_number)
    if some_number >=0 && some_number <= 1
        sf = some_number;
    else
        error(['Scale factor ' num2str(some_number) ' should be positive.'])
    end
end
end