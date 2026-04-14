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
% --------------------------------------------------------------------------

%% muscle
% Find out which muscles span which joint, thus interacts with its associated coordinates.
model_info.muscle_info.muscle_spanning_joint_info = get_muscle_spanning_joint_info(S,osim_path,model_info);

% Coordinates actuated by muscles
model_info.ExtFunIO.jointi.muscleActuated = find(sum(model_info.muscle_info.muscle_spanning_joint_info,1)>0);


% Use polynomial approximatrion
if strcmp(S.misc.msk_geom_eq,'polynomials') 
    % Only perform muscle analysis and fitting if the result is not yet
    % available, because the analysis takes long. (4 minutes)

    msk_geom_path = fullfile(S.misc.subject_path,[S.misc.msk_geom_name,'.casadi']);

    if ~isfile(msk_geom_path) || S.misc.msk_geom_always_new_fit
        % Analyze the muscle-tendon lengths, velocities, and moment arms in 
        % function of coordinate values
        t0 = tic;
        disp(['   start analysing musculoskeletal geometry...'])
        muscle_data = muscleAnalysisAPI(S,osim_path,model_info);
        disp(['   analysing MSK geometry: ' num2str(toc(t0),'%.2f') ' s'])

        % fit polynomial to approximate the results
        t1 = tic;
        model_info.muscle_info.polyFit.MuscleInfo = PolynomialFit(S,...
            muscle_data, model_info.muscle_info.muscle_spanning_joint_info);
        disp(['   approximating MSK geometry: ' num2str(toc(t1),'%.2f') ' s'])
%         disp(['   (total duration: ' num2str(toc(t0)) ' s)'])
        
        % Save sampling and fitting data
        msk_geom_fit_info.samples = muscle_data;
        msk_geom_fit_info.fit = model_info.muscle_info.polyFit.MuscleInfo;
        save(fullfile(S.misc.subject_path,[S.misc.msk_geom_name,'_info.mat']),...
            'msk_geom_fit_info')
        
    else
        disp(['   using existing musculoskeletal geometry ' S.misc.msk_geom_name])
    end

else
    % Other fitting methods are not (yet) implemented
    warning(['Selected method to approximate musculoskeletal geometry: "',...
        S.misc.msk_geom_eq, '" not implemented.'])
end

%% ligaments
if model_info.ligament_info.NLigament > 0
    % Find out which ligaments span which joint, thus interacts with its 
    % associated coordinates.
    model_info.ligament_info.ligament_spanning_joint_info = ...
        get_muscle_spanning_joint_info(S,osim_path,model_info,'ligaments');
    
    % Separate out ligaments that only cross 1 coordinate
    n_coordinates_crossed = sum(model_info.ligament_info.ligament_spanning_joint_info,2);
    model_info.ligament_info.ligament_spanning_single_coord = ...
        model_info.ligament_info.ligament_spanning_joint_info;
    model_info.ligament_info.ligament_spanning_single_coord(n_coordinates_crossed(:)>1,:) = 0;
    model_info.ligament_info.ligament_spanning_multi_coord = ...
        model_info.ligament_info.ligament_spanning_joint_info;
    model_info.ligament_info.ligament_spanning_multi_coord(n_coordinates_crossed(:)==1,:) = 0;
    
    
    
    % Use polynomial approximation
    if strcmp(S.misc.msk_geom_eq,'polynomials') 
        
        % Analyze the ligament lengths, velocities, and moment arms in 
        % function of coordinate values
        t0 = tic;
        ligament_data = muscleAnalysisAPI(S,osim_path,model_info,'ligaments');
        disp(['   analysing ligament geometry: ' num2str(toc(t0),'%.2f') ' s'])
    
        model_info.ligament_info.polyFit.DummySamples = ligament_data;

        % fit polynomial to approximate the results
        if sum(model_info.ligament_info.ligament_spanning_multi_coord,'all') > 0
            t1 = tic;
            model_info.ligament_info.polyFit.LigamentInfo = PolynomialFit(S,...
                ligament_data, model_info.ligament_info.ligament_spanning_multi_coord);
            
            disp(['   approximating ligament geometry: ' num2str(toc(t1),'%.2f') ' s'])
        end
    
    else
        % Other fitting methods are not (yet) implemented
        warning(['Selected method to approximate musculoskeletal geometry: "',...
            S.misc.msk_geom_eq, '" not implemented.'])
    end

end



end % end of function
