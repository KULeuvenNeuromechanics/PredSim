function [MuscleInfo] = PolynomialFit(S,MuscleData,muscle_spanning_joint_info)
% --------------------------------------------------------------------------
% PolynomialFit
%   This function computes the polynomials to approximate muscle-tendon 
%   lengths, velocities and moment arms. 
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - MuscleData -
%   * structure with joint angles and according muscle-tendon lengths and
%       momentarms (output from muscleAnalysisAPI.m)
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - MuscleInfo -
%   * structure with information about the polynomials
% 
% Original author: Lars D'Hondt
% Original date: 15 October 2025
% --------------------------------------------------------------------------


%%
% On HPC (any linux?):
% When using OpenBLAS as BLAS/LAPACK libraries, the built-in mldivide
% leads to segmentation violations. The current workaround is to provide
% a custom function in the LinearAlgebra subdirectory that directly calls
% LAPACK, which is compiled below if not yet done.
lapack_version = version('-lapack');
if startsWith(lapack_version, 'Intel')
    % Use built-in matlab function
    mldivide_impl = 'mldivide';
else
    % Compile if not yet done
    if (exist('mldivide_lapack') == 0)
        currentFolder = pwd;
        cd(fullfile(S.misc.main_path, 'LinearAlgebra'))
        mex mldivide_lapack.c -lopenblas
        cd(currentFolder)
    end
    mldivide_impl = 'mldivide_lapack';
end


%% Construct the polynomials for the moment arms and muscle length

for m_nr=1:size(MuscleData.lMT,2)
    
    index_dof_crossing = find(muscle_spanning_joint_info(m_nr,:)==1);

    if isempty(index_dof_crossing)
        continue
    end
    
    [coeff, stats, mu] = mvpolyfit(MuscleData.q(:,index_dof_crossing),...
        MuscleData.lMT(:,m_nr),...
        [S.misc.poly_order.lower, S.misc.poly_order.upper],...
        -squeeze(MuscleData.dM(:,m_nr,index_dof_crossing)),... % dM = -jac(lMT,q)
        "reduced_coeff",S.misc.reduce_coeff_fit,...
        "threshold_rmse_y",S.misc.threshold_lMT_fit,...
        "threshold_rmse_ydx",S.misc.threshold_dM_fit,...
        "mldivide_impl",mldivide_impl);
    
    
    MuscleInfo.muscle(m_nr).DOF = MuscleData.dof_names(index_dof_crossing);
    MuscleInfo.muscle(m_nr).m_name = MuscleData.muscle_names{m_nr};
    MuscleInfo.muscle(m_nr).coeff = coeff;
    MuscleInfo.muscle(m_nr).order = stats.order;
    MuscleInfo.muscle(m_nr).lMT_error_rms = stats.rmse_y;
    MuscleInfo.muscle(m_nr).dm_error_rms = stats.rmse_ydx;
    MuscleInfo.muscle(m_nr).stats = stats;
    MuscleInfo.muscle(m_nr).mu = mu;
                
end

% figure();
% tiledlayout('flow')
% nexttile
% hold on;
% plot(lMT_all_error)
% xlimits = get(gca, 'XLim');
% plot(xlimits, [1, 1]*S.misc.threshold_lMT_fit, 'r', 'linewidth', 2)
% title('RMS error on the approximated muscle-tendon length')
% ylabel('RMS error (m)')
% 
% nexttile
% hold on;
% plot(max(DM_all_error, [], 2))
% xlimits = get(gca, 'XLim');
% plot(xlimits, [1, 1]*S.misc.threshold_dM_fit, 'r', 'linewidth', 2)
% title('maximal RMS error on the approximated muscle moment arm')
% ylabel('RMS error (m)')
% 
% nexttile
% hold on;
% plot(order_all)
% ylim([0 max_order+1])
% xlimits = get(gca, 'XLim');
% plot(xlimits, [max_order, max_order], 'r', 'linewidth', 2)
% title('Order of the polynomial approximation')
% ylabel('Order')


end

