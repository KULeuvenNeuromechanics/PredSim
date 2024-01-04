function [f_ligamentMoment,f_ligamentMoment_single,f_ligamentMoment_multi,...
    f_ligamentLengthForce] = createCasadi_Ligaments(S,model_info)
% --------------------------------------------------------------------------
% createCasadi_ContractDynam
%   Function to create Casadi functions for muscle contraction dynamics.
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

Qs_MX = MX.sym('Qs',model_info.ExtFunIO.jointi.nq.all,1);
Qdots_MX = MX.sym('Qdots',model_info.ExtFunIO.jointi.nq.all,1);

% initialise moments from ligaments with vectors of sparse zeros
% ligaments that cross a single degree of freedon
M_lig_single_MX = MX(model_info.ExtFunIO.jointi.nq.all,1);
% ligaments that cross multiple degrees of freedom
M_lig_multi_MX = MX(model_info.ExtFunIO.jointi.nq.all,1);

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
            % test for valid input force-length
            if ~exist(model_info.ligament_info.parameters(j).stiffness,'file')
                error(['Unable to find ' model_info.ligament_info.parameters(j).stiffness,...
                    '. Add the function to ./ModelComponents, or select an existing function for ',...
                    'the force-length of ligament ' model_info.ligament_info.parameters(j).ligament_name '.'])
            end
    
            % evaluate force-length for dummy motion
            f_stiffness = str2func(model_info.ligament_info.parameters(j).stiffness);
            CSA = model_info.ligament_info.parameters(j).cross_section_area;
            ls = model_info.ligament_info.parameters(j).slack_length;
            l_lig = model_info.ligament_info.polyFit.DummySamples.lMT(idx_sort,j);
            d_lig = model_info.ligament_info.polyFit.DummySamples.dM(idx_sort,j,i);
            F_lig_j = f_stiffness(CSA,ls,l_lig);
    
            % ligament length and force for use in post-processing
            f_L_lig_j = interpolant(['f_L_lig_' model_info.ligament_info.parameters(j).ligament_name],...
                'bspline',{q_i},l_lig);
            f_F_lig_j = interpolant(['f_F_lig_' model_info.ligament_info.parameters(j).ligament_name],...
                'bspline',{q_i},F_lig_j);
            
            L_lig_pp_MX(j) = f_L_lig_j(Qs_MX(i));
            F_lig_pp_MX(j) = f_F_lig_j(Qs_MX(i));


            % calculate moment
            M_lig_j = F_lig_j.*d_lig;
            M_i = M_i + M_lig_j;
        end
    
        % create look-up table
        f_lig_single_coord = interpolant(['f_lig_single_coord_',...
            model_info.ExtFunIO.coord_names.all{i}], 'bspline',{q_i},M_i);
        
        % assign to coordinate
        M_lig_single_MX(i) = f_lig_single_coord(Qs_MX(i));
    
    end
    
    
    %% Ligaments that cross multiple coordinates
    idx_multi_coord = find(sum(model_info.ligament_info.ligament_spanning_multi_coord,2));

    % Assemble polynomial approximation from coefficients
    qin_SX = SX.sym('qin',model_info.ExtFunIO.jointi.nq.all,1);
    l_lig_p_SX = SX(model_info.ligament_info.NLigament,1);

    if sum(model_info.ligament_info.ligament_spanning_multi_coord,'all') > 0
        ligament_spanning_info_m = model_info.ligament_info.ligament_spanning_multi_coord(:,:);
        LigamentInfo_m.muscle = model_info.ligament_info.polyFit.LigamentInfo.muscle(:);
        load(fullfile(S.misc.main_path,'CasadiFunctions','nCoeffMat.mat'), 'nCoeffMat');
        load(fullfile(S.misc.main_path,'CasadiFunctions','expoVal_all.mat'), 'expoVal_all');
        for i=1:model_info.ligament_info.NLigament      
            index_dof_crossing = find(ligament_spanning_info_m(i,:)==1);
            nr_dof_crossing = length(index_dof_crossing);
            if nr_dof_crossing == 0
                continue
            end
            order = LigamentInfo_m.muscle(i).order;
            % get matrix with variabes for polynomial
            % (e.g. mat = [1, q_1, q_1^2] for a 2nd order polynomial in 1 variable)
            [mat,~] = n_art_mat_3_cas_SX_7(qin_SX(index_dof_crossing,1)',...
                nCoeffMat(order,nr_dof_crossing), expoVal_all{order,nr_dof_crossing});
            % multiply variable matrix with coefficient vactor to get polynomial expression
            l_lig_p_SX(i,1) = mat'*LigamentInfo_m.muscle(i).coeff; 
        end 
    end
    f_l_lig = Function('f_l_lig',{qin_SX},{l_lig_p_SX});
    
    % Create casadi function force-length
    l_lig_SX = SX.sym('l_lig',model_info.ligament_info.NLigament,1);
    F_lig_SX = SX(model_info.ligament_info.NLigament,1);
    for j=idx_multi_coord'
        % test for valid input force-length
        if ~exist(model_info.ligament_info.parameters(j).stiffness,'file')
            error(['Unable to find ' model_info.ligament_info.parameters(j).stiffness,...
                '. Add the function to ./ModelComponents, or select an existing function for ',...
                'the force-length of ligament ' model_info.ligament_info.parameters(j).ligament_name '.'])
        end
        % evaluate force-length
        f_stiffness = str2func(model_info.ligament_info.parameters(j).stiffness);
        CSA = model_info.ligament_info.parameters(j).cross_section_area;
        ls = model_info.ligament_info.parameters(j).slack_length;
        F_lig_SX(j) = f_stiffness(CSA,ls,l_lig_SX(j));
    end
    
    f_F_lig = Function('f_F_lig',{l_lig_SX},{F_lig_SX});
    
    % Calculate moments
    l_lig_MX = f_l_lig(Qs_MX);
    F_lig_MX = f_F_lig(l_lig_MX);
    % transposed Jacobian of lengths to angles (i.e. moment arms), multiplied by forces
    M_lig_multi_MX = -jtimes(l_lig_MX,Qs_MX,F_lig_MX,true); 

    % Lengths and forces for post-processing
    L_lig_pp_MX = L_lig_pp_MX + l_lig_MX; % both vectors cannot have a non-zero on the same index.
    F_lig_pp_MX = F_lig_pp_MX + F_lig_MX; 
    
    v_lig_pp_MX = jtimes(L_lig_pp_MX,Qs_MX,Qdots_MX);

else
    L_lig_pp_MX = MX(1,1);
    v_lig_pp_MX = MX(1,1);
    F_lig_pp_MX = MX(1,1);
end

% total
M_lig_MX = M_lig_single_MX + M_lig_multi_MX;

f_ligamentMoment = Function('f_ligamentMoment',{Qs_MX},{M_lig_MX},{'Qs'},{'M_lig'});
f_ligamentMoment_single = Function('f_ligamentMoment_single',{Qs_MX},{M_lig_single_MX},...
    {'Qs'},{'M_lig'});
f_ligamentMoment_multi = Function('f_ligamentMoment_multi',{Qs_MX},{M_lig_multi_MX},...
    {'Qs'},{'M_lig'});
f_ligamentLengthForce = Function('f_ligamentLengthForce',{Qs_MX,Qdots_MX},{L_lig_pp_MX,...
    v_lig_pp_MX F_lig_pp_MX},{'Qs','Qdots'},{'length','velocity','force'});


end % end of function