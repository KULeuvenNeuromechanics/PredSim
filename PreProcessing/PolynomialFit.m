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
% Original author: Original code from Wouter Aerts, adapted by Antoine Falisse
% Original date: 19/December/2018
%
% update:
%   Compatibility with generalized code structure. 
%   Get muscle_spanning_joint_info from input 
%
% Last edit by: Lars D'Hondt
% Last edit date: 05/April/2022
% --------------------------------------------------------------------------


%%
% When using OpenBLAS as BLAS/LAPACK libraries, the built-in mldivide
% leads to segmentation violations. The current workaround is to provide
% a custom function in the LinearAlgebra subdirectory that directly calls
% LAPACK, which is compiled below if not yet done.
lapack_version = version('-lapack');
if startsWith(lapack_version, 'Intel')
    mldivide_impl = [];
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

muscle_sel = 1:length(MuscleData.muscle_names);

q_all = MuscleData.q;

max_order = S.misc.poly_order.upper;
nr_samples = length(q_all(:,1));

lMT_all_error = zeros(length(muscle_sel), 1);
DM_all_error = zeros(length(muscle_sel), length(q_all(1,:)));
order_all = zeros(length(muscle_sel), 1);

for m_nr=1:length(muscle_sel)
    muscle_index = muscle_sel(m_nr);
    
    index_dof_crossing = find(muscle_spanning_joint_info(muscle_index,:)==1);
    nr_dof_crossing = length(index_dof_crossing);
    if nr_dof_crossing == 0
        continue
    end
    
    lMT = MuscleData.lMT(:,muscle_index);
    dM = zeros(nr_samples, nr_dof_crossing);
    for dof_nr = 1:nr_dof_crossing
        dM(:,dof_nr) = MuscleData.dM(:,muscle_index,index_dof_crossing(dof_nr));
    end

    [coeff, stats, mu] = mvpolyfit(q_all(:,index_dof_crossing), lMT,...
        [S.misc.poly_order.lower,S.misc.poly_order.upper], -dM,...
        "reduced_coeff",1,...
        "threshold_rmse_y",S.misc.threshold_lMT_fit,...
        "threshold_rmse_ydx",S.misc.threshold_dM_fit,...
        "mldivide",mldivide_impl);
    
    
    MuscleInfo.muscle(m_nr).DOF = MuscleData.dof_names(index_dof_crossing);
    MuscleInfo.muscle(m_nr).m_name = MuscleData.muscle_names{muscle_index};
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

