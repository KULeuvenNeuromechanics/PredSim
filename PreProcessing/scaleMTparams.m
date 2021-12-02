function [model_info] = scaleMTparams(S,model_info)
% --------------------------------------------------------------------------
%scaleMTparams 
%   scale the MT parameters based on the input of the user in the settings
%   and store the scaled values in model_info.
% 
% INPUT:
%   - S -
%   * Settings structure, including scaling settings for the muscle tendon
%   parameters.
%  
%   - model_info -
%   * The info of the model read from the OpenSim model.
% 
% OUTPUT:
%   - model_info -
%   * Update model_info struct where the old MT parameter values are
%   replaced by the scaled values.
% 
% Original author: Bram Van Den Bosch
% Original date: 02/12/2021
%
% Last edit by: Bram Van Den Bosch
% Last edit date: 02/12/2021
% --------------------------------------------------------------------------

%% read MT params from model_info
Fmax      = model_info.muscle_info.params.Fmax;
% lMopt     = model_info.muscle_info.params.lMopt;
% lTs       = model_info.muscle_info.params.lTs;
% pennAngle = model_info.muscle_info.params.pennAngle;
% Vmax      = model_info.muscle_info.params.Vmax;

muscle_names = model_info.muscle_info.params.names;

%% read scaling paramaters from S
Fmax_scale_names   = S.subject.muscle_strength.names;
Fmax_scale_factors = S.subject.muscle_strength.scale_factors;
% lMopt_scale     =
% lTs_scale       = 
% pennAngle_scale = 
% Vmax_scale      = 


%% perform scaling of MT params
for i = 1:length(Fmax_scale_names)
    idx = find(strcmp(Fmax_scale_names(i),muscle_names(:,:)));
    Fmax(idx) = Fmax(idx)*Fmax_scale_factors(i);
end

%% write new 
% model_info.muscle_info.params.Fmax = Fmax;

end

