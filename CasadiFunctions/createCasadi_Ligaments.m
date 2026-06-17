function [f_ligamentMoment,f_ligamentMoment_single,f_ligamentMoment_multi,...
    f_ligamentLengthForce] = createCasadi_Ligaments(S,model_info)
% --------------------------------------------------------------------------
% createCasadi_ContractDynam
%   Function to create Casadi functions for ligament dynamics.
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - f_ligamentMoment -
%   * function for angle-moment of all ligaments
%
%   - f_ligamentMoment_single -
%   * function for angle-moment of ligaments that cross 1 coordinate
%
%   - f_ligamentMoment_multi -
%   * function for angle-moment of ligaments that cross multiple
%   coordinates
% 
%   - f_ligamentLengthForce -
%   * function for length, lenthening speed, and force of ligaments in function of joint
%   positions
%
%
% Original author: Lars D'Hondt
% Original date: 5/April/2023
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.*

if isfile(fullfile(S.misc.subject_path,S.misc.lig_force_length_name)) && isfile(fullfile(S.misc.subject_path,S.misc.lig_momTot_name))

    f_ligamentMoment = Function.load(fullfile(S.misc.subject_path,S.misc.lig_momTot_name));
    f_ligamentMoment_single = Function.load(fullfile(S.misc.subject_path,S.misc.lig_momSin_name));
    f_ligamentMoment_multi = Function.load(fullfile(S.misc.subject_path,S.misc.lig_momMulti_name));
    f_ligamentLengthForce = Function.load(fullfile(S.misc.subject_path,S.misc.lig_force_length_name));

else

    Qs_MX = MX.sym('Qs',model_info.ExtFunIO.jointi.nq.all,1);
    Qdots_MX = MX.sym('Qdots',model_info.ExtFunIO.jointi.nq.all,1);
    
    M_lig_single_MX = MX(model_info.ExtFunIO.jointi.nq.all,1);
    
    if model_info.ligament_info.NLigament > 0
    
        L_lig_pp_MX = MX(model_info.ligament_info.NLigament,1);
        F_lig_pp_MX = MX(model_info.ligament_info.NLigament,1);
    
        %% Ligaments that cross a single coordinate
        % These are grouped as a lumped angle-torque
    
        idx_coord_single = find(sum(model_info.ligament_info.ligament_spanning_single_coord,1));
        
        for i=idx_coord_single
            q_i = model_info.ligament_info.polyFit.DummySamples.q(:,i);
            M_i = zeros(size(q_i));
        
            idx_lig = find(model_info.ligament_info.ligament_spanning_single_coord(:,i));
        
            % sort angle and moment in order of angle
            [q_i,idx_sort] = sort(q_i);
    
            for j=idx_lig'
       
                % evaluate force-length for dummy motion
                try
                    f_stiffness = str2func(model_info.ligament_info.parameters(j).stiffness);
                catch
                    f_stiffness = str2func(char(model_info.ligament_info.parameters(j).stiffness));
                end
                
                linStiff = model_info.ligament_info.parameters(j).linear_stiffness;
                ls = model_info.ligament_info.parameters(j).slack_length;
                l_lig = model_info.ligament_info.polyFit.DummySamples.lMT(idx_sort,j);
                d_lig = model_info.ligament_info.polyFit.DummySamples.dM(idx_sort,j,i);
                F_lig_j = f_stiffness(linStiff,ls,l_lig)';
        
                % ligament length and force for use in post-processing
                f_L_lig_j = interpolant(['f_L_lig_' char(model_info.ligament_info.parameters(j).ligament_name)],...
                    'bspline',{q_i},l_lig);
                f_F_lig_j = interpolant(['f_F_lig_' char(model_info.ligament_info.parameters(j).ligament_name)],...
                    'bspline',{q_i},F_lig_j);
                
                L_lig_pp_MX(j) = f_L_lig_j(Qs_MX(i));
                F_lig_pp_MX(j) = f_F_lig_j(Qs_MX(i));
    
    
                % calculate moment
                M_lig_j = F_lig_j.*d_lig;
                M_i = M_i + M_lig_j;
            end
        
            % create look-up table
            f_lig_single_coord = interpolant(['f_lig_single_coord_',...
                char(model_info.ExtFunIO.coord_names.all{i})], 'bspline',{q_i},M_i);
            
            % assign to coordinate
            M_lig_single_MX(i) = f_lig_single_coord(Qs_MX(i));
        
        end
        
        
        %% Ligaments that cross multiple coordinates
        
        %% Lengths 
        idx_multi_coord = find(sum(model_info.ligament_info.ligament_spanning_multi_coord,2));
    
        % Assemble polynomial approximation from coefficients
        qin_SX = SX.sym('qin',model_info.ExtFunIO.jointi.nq.all,1);
        l_lig_p_SX = SX(model_info.ligament_info.NLigament,1);
    
        if sum(model_info.ligament_info.ligament_spanning_multi_coord,'all') > 0
            ligament_spanning_info_m = model_info.ligament_info.ligament_spanning_multi_coord(:,:);
            LigamentInfo_m.muscle = model_info.ligament_info.polyFit.LigamentInfo.muscle(:);
            load nCoeffMat
            load expoVal_all
            for i=1:model_info.ligament_info.NLigament      
                index_dof_crossing = find(ligament_spanning_info_m(i,:)==1);
                nr_dof_crossing = length(index_dof_crossing);
                if nr_dof_crossing == 0
                    continue
                end
                order = LigamentInfo_m.muscle(i).order;
                % original 
                %[mat,~] = n_art_mat_3_cas_SX_7(qin_SX(index_dof_crossing,1)',...
                %    nCoeffMat(order,nr_dof_crossing), expoVal_all{order,nr_dof_crossing});
                %Gil imp
                [mat,~]    = n_art_mat_9_GC_SX(qin_SX(index_dof_crossing,1)',...
                     order);

                mat = mat';

                l_lig_p_SX(i,1) = mat'*LigamentInfo_m.muscle(i).coeff; 
            end 
        end
         f_l_lig = Function('f_l_lig',{qin_SX},{l_lig_p_SX});
        
        
        %% Forces 
        idx_multi_coord = find(sum(model_info.ligament_info.ligament_spanning_multi_coord,2));
    
        % Assemble polynomial approximation from coefficients
        qin_SX = SX.sym('qin',model_info.ExtFunIO.jointi.nq.all,1);
        f_lig_p_SX = SX(model_info.ligament_info.NLigament,1);
    
        if sum(model_info.ligament_info.ligament_spanning_multi_coord,'all') > 0
            ligament_spanning_info_m = model_info.ligament_info.ligament_spanning_multi_coord(:,:);
            LigamentInfo_m.muscle = model_info.ligament_info.polyFit.LigamentForceInfo.muscle(:);
            load nCoeffMat % I have no idea what this measn
            load expoVal_all % what even are these ?
            for i=1:model_info.ligament_info.NLigament      
                index_dof_crossing = find(ligament_spanning_info_m(i,:)==1);
                nr_dof_crossing = length(index_dof_crossing);
                if nr_dof_crossing == 0
                    continue
                end
                order = LigamentInfo_m.muscle(i).order;
                % original
                %[mat,~] = n_art_mat_3_cas_SX_7(qin_SX(index_dof_crossing,1)',...
                %    nCoeffMat(order,nr_dof_crossing), expoVal_all{order,nr_dof_crossing});

                % Gil 
                [mat,~] = n_art_mat_9_GC_SX(qin_SX(index_dof_crossing,1)',...
                    order);
                mat = mat'; 

                f_lig_p_SX(i,1) = mat'*LigamentInfo_m.muscle(i).coeff; 
            end 
        end
         f_F_lig = Function('f_F_lig',{qin_SX},{f_lig_p_SX});
        
        
         %%
    
        % Calculate moments
        l_lig_MX = f_l_lig(Qs_MX);
        % If I change the input here to be Qs_MX this means I get force from Qs
        % right ? 
        %F_lig_MX = f_F_lig(l_lig_MX);
        F_lig_MX = f_F_lig(Qs_MX);
        % transposed Jacobian of lengths to angles (i.e. moment arms), multiplied by forces
        M_lig_multi_MX = -jtimes(l_lig_MX,Qs_MX,F_lig_MX,true); 
    
        % Lengths and forces for post-processing
        L_lig_pp_MX = L_lig_pp_MX + l_lig_MX; % both vectors cannot have a non-zero on the same index.
        F_lig_pp_MX = F_lig_pp_MX + F_lig_MX; 
        
        v_lig_pp_MX = jtimes(L_lig_pp_MX,Qs_MX,Qdots_MX);
    
    else
        % if there are no ligaments, moments are zero (sparse)
        M_lig_multi_MX = MX(model_info.ExtFunIO.jointi.nq.all,1);
    end
    
    % total
    M_lig_MX = M_lig_single_MX + M_lig_multi_MX;
    

    f_ligamentMoment = Function('f_ligamentMoment',{Qs_MX},{M_lig_MX},{'Qs'},{'M_lig'});
    f_ligamentMoment.save(fullfile(S.misc.subject_path,S.misc.lig_momTot_name));

    f_ligamentMoment_single = Function('f_ligamentMoment_single',{Qs_MX},{M_lig_single_MX},...
        {'Qs'},{'M_lig'});
    f_ligamentMoment_single.save(fullfile(S.misc.subject_path,S.misc.lig_momSin_name));

    f_ligamentMoment_multi = Function('f_ligamentMoment_multi',{Qs_MX},{M_lig_multi_MX},...
        {'Qs'},{'M_lig'});
    f_ligamentMoment_multi.save(fullfile(S.misc.subject_path,S.misc.lig_momMulti_name));

    f_ligamentLengthForce = Function('f_ligamentLengtForce',{Qs_MX,Qdots_MX},{L_lig_pp_MX,...
        v_lig_pp_MX F_lig_pp_MX},{'Qs','Qdots'},{'length','velocity','force'});
    f_ligamentLengthForce.save(fullfile(S.misc.subject_path,S.misc.lig_force_length_name));
    
end

end % end of function