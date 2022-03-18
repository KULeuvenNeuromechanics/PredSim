function [model_info] = get_musculoskeletal_geometry_approximation(S,osim_path,model_info)
% --------------------------------------------------------------------------
% get_musculoskeletal_geometry_approximation
%   Analyzes the muscle-tendon lengths, velocities, and moment arms in
%   function of coordinate values, and fits expressions to approximate the
%   results.
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
% Original date: 18/March/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

% Analyze the muscle-tendon lengths, velocities, and moment arms in function of coordinate values
muscle_data = muscleAnalysis(S,osim_path,model_info);

% fit expressions to approximate the results
if strcmp( S.misc.msk_geom_eq,'polynomials')
    [model_info] = PolynomialFit(S,muscle_data);
else
    [model_info] = PolynomialFit(S,muscle_data);
    warning(['Selected method to approximate musculoskeletal geometry not',...
        ' implemented, using polynomials instead.'])
end
