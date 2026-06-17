function [MuscleInfo] = PolynomialFit_BlankevoortLigament(S,MuscleData,muscle_spanning_joint_info, model_info)
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
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% Original author: Original code from Wouter Aerts, adapted by Antoine Falisse
% Original date: 19/December/2018
%
% update:
%   Compatibility with generalized code structure. 
%   Get muscle_spanning_joint_info from input 
%
% Last edit by: Bryce Killen
% Last edit date: 08/August/2023
% --------------------------------------------------------------------------


muscle_sel=[];

for m_nr = 1:length(MuscleData.muscle_names)
    if strcmp(MuscleData.muscle_names{m_nr}(end-1:end), '_r')...
            || strcmp(MuscleData.muscle_names{m_nr}(end-1:end), '_l') % was _l before
        muscle_sel = [muscle_sel m_nr];
    end
end

%muscle_spanning_joint_info = model_info.muscle_info.muscle_spanning_joint_info;
q_all = MuscleData.q;

max_order = S.misc.poly_order.upper;
nr_samples = length(q_all(:,1));

fMT_all_error = zeros(length(muscle_sel), 1);
order_all = zeros(length(muscle_sel), 1);

for m_nr=1:length(muscle_sel)
    muscle_index = muscle_sel(m_nr);
    
    index_dof_crossing = find(muscle_spanning_joint_info(muscle_index,:)==1);
    nr_dof_crossing = length(index_dof_crossing);
    
    if nr_dof_crossing == 0
        continue
    end
    
    lMT = MuscleData.lMT(:,muscle_index);
    % For the set of lMT - get a set of forces 
    % get properties of the ligaments 
    linStiff = model_info.ligament_info.parameters(m_nr).linear_stiffness;
    ls = model_info.ligament_info.parameters(m_nr).slack_length;
    fMT = BlankevoortLig(linStiff,ls,lMT)';
       
    criterion_full_filled = 0;
    order = S.misc.poly_order.lower;
    while criterion_full_filled==0
        % TRY TO SUB IN THE CODE OF GIL
        %ORIG
        %[mat,diff_mat_q] = n_art_mat_3(q_all(:,index_dof_crossing), order);
        % GIL - after removing the ITB - in theory the old should work
        [mat,diff_mat_q] = n_art_mat_9_GC(q_all(:,index_dof_crossing), order);

        nr_coeffs = length(mat(1,:));
        
        %diff_mat_q_all = zeros(nr_samples*nr_dof_crossing, nr_coeffs);
        %for dof_nr = 1:nr_dof_crossing
        %    diff_mat_q_all(nr_samples*(dof_nr-1)+1:nr_samples*dof_nr,:) =...
        %        -squeeze(diff_mat_q(:,:,dof_nr));
        %end
        
        %coeff=[mat ; diff_mat_q_all]\[lMT; dMT(:)];
        coeff=[mat]\[fMT];
        %dM_recon = zeros(nr_samples, nr_dof_crossing);
        %for dof_nr = 1:nr_dof_crossing
        %    dM_recon(:,dof_nr) = (-squeeze(diff_mat_q(:,:,dof_nr)))*coeff;
        %end
        %lMT_recon=mat*coeff;
        fMT_recon=mat*coeff;
        
        fMT_error_rms = sqrt(mean((fMT - fMT_recon).^2));
        %lMT_error_rms = sqrt(mean((lMT - lMT_recon).^2));
        %dm_error_rms = sqrt(mean((dM - dM_recon).^2));
        
        %criterion_full_filled = lMT_error_rms<=S.misc.threshold_lMT_fit ...
        %    & max(dm_error_rms)<=S.misc.threshold_dM_fit;
        
        criterion_full_filled = fMT_error_rms<=1;
        
        if order==max_order
            criterion_full_filled = 1;
        end
        if criterion_full_filled==0
            order = order+1;
        end
    end
    
    MuscleInfo.muscle(m_nr).DOF = MuscleData.dof_names(index_dof_crossing);
    MuscleInfo.muscle(m_nr).m_name = MuscleData.muscle_names{muscle_index};
    MuscleInfo.muscle(m_nr).coeff = coeff;
    MuscleInfo.muscle(m_nr).order = order;
    %MuscleInfo.muscle(m_nr).lMT_error_rms = lMT_error_rms;
    %MuscleInfo.muscle(m_nr).dm_error_rms = dm_error_rms;
    
    MuscleInfo.muscle(m_nr).fMT_error_rms = fMT_error_rms;
    
    %lMT_all_error(m_nr) = lMT_error_rms;
    fMT_all_error(m_nr) = fMT_error_rms;
    %DM_all_error(m_nr, index_dof_crossing) = dm_error_rms;
    
    order_all(m_nr) = order;            
end


end

