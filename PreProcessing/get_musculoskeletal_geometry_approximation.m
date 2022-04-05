function [model_info] = get_musculoskeletal_geometry_approximation(S,...
    osim_path,model_info)
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

% Find out which muscles span wich joint, thus interacts with its associated coordinates.
model_info.muscle_info.muscle_spanning_joint_info = get_muscle_spanning_joint_info(S,osim_path,model_info);

if ~strcmp(S.misc.msk_geom_eq,'polynomials') 
    warning(['Selected method to approximate musculoskeletal geometry: "',...
        S.misc.msk_geom_eq, '" not implemented, using polynomials instead.'])
    S.misc.msk_geom_eq = 'polynomials';
end

% Use polynomial approximatrion
if strcmp(S.misc.msk_geom_eq,'polynomials') 
    % Only perform muscle analysis and fitting if the result is not yet
    % available, because the analysis takes long. (4 minutes)
    if ~isfile(fullfile(S.misc.subject_path,'f_lMT_vMT_dM_poly'))
        % Analyze the muscle-tendon lengths, velocities, and moment arms in function of coordinate values
        t0 = tic;
        % muscle_data = muscleAnalysis(S,osim_path,model_info); % to be removed since too slow
        muscle_data = muscleAnalysisAPI(S,osim_path,model_info); % faster version
        disp(['   analysing MSK geometry: ' num2str(toc(t0)) ' s'])

        % fit polynomial to approximate the results
        t1 = tic;
        [model_info] = PolynomialFit(S,muscle_data,model_info);
        disp(['   approximating MSK geometry: ' num2str(toc(t1)) ' s'])
%         disp(['total duration: ' num2str(toc(t0)) ' s'])
    end

else
    % Other fitting methods are not (yet) implemented
end



