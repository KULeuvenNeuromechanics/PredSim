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

% Find out which muscles span which joint, thus interacts with its associated coordinates.
model_info.muscle_info.muscle_spanning_joint_info = get_muscle_spanning_joint_info(S,osim_path,model_info);

% Coordinates actuated by muscles
model_info.ExtFunIO.jointi.muscleActuated = find(sum(model_info.muscle_info.muscle_spanning_joint_info,1)>0);


% Use polynomial approximatrion
if strcmp(S.misc.msk_geom_eq,'polynomials') 
    % Only perform muscle analysis and fitting if the result is not yet
    % available, because the analysis takes long. (4 minutes)
    if ~isfile(fullfile(S.misc.subject_path,S.misc.msk_geom_name)) || isempty(S.misc.msk_geom_bounds)
        % Analyze the muscle-tendon lengths, velocities, and moment arms in function of coordinate values
        t0 = tic;
        muscle_data = muscleAnalysisAPI(S,osim_path,model_info); % faster version
        disp(['   analysing MSK geometry: ' num2str(toc(t0)) ' s'])

        % fit polynomial to approximate the results
        t1 = tic;
        [model_info] = PolynomialFit(S,muscle_data,model_info);
        disp(['   approximating MSK geometry: ' num2str(toc(t1)) ' s'])
        disp(['   (total duration: ' num2str(toc(t0)) ' s)'])
    end

else
    % Other fitting methods are not (yet) implemented
    warning(['Selected method to approximate musculoskeletal geometry: "',...
        S.misc.msk_geom_eq, '" not implemented.'])
end



