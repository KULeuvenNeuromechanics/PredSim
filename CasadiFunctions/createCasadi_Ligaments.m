function [f_ligamentMoment,f_ligamentMoment_single,f_ligamentMoment_multi] = createCasadi_Ligaments(S,model_info)
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
% Original author: Lars D'Hondt
% Original date: 5/April/2023
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.*

Qs_MX = MX.sym('Qs',model_info.ExtFunIO.jointi.nq.all,1);
M_lig_single_MX = MX(model_info.ExtFunIO.jointi.nq.all,1);

if model_info.ligament_info.NLigament > 0

    %% Ligaments that cross a single coordinate
    % These are grouped as a lumped angle-torque

    idx_coord_single = find(sum(model_info.ligament_info.ligament_spanning_single_coord,1));
    
    for i=idx_coord_single
        q_i = model_info.ligament_info.polyFit.DummySamples.q(:,i);
        M_i = zeros(size(q_i));
    
        idx_lig = find(model_info.ligament_info.ligament_spanning_single_coord(:,i));
    
        for j=idx_lig'
            % test for valid input force-length
            if ~exist(model_info.ligament_info.parameters(j).stiffness,'file')
                error(['Unable to find ' model_info.ligament_info.parameters(j).stiffness,...
                    '. Add the function to ./ModelComponents, or select an existing function for ',...
                    'the force-length of ligament ' model_info.ligament_info.parameters(j).name '.'])
            end
    
            % evaluate force-length for dummy motion
            f_stiffness = str2func(model_info.ligament_info.parameters(j).stiffness);
            CSA = model_info.ligament_info.parameters(j).cross_section_area;
            ls = model_info.ligament_info.parameters(j).slack_length;
            l_lig = model_info.ligament_info.polyFit.DummySamples.lMT(:,j);
            d_lig = model_info.ligament_info.polyFit.DummySamples.dM(:,j,i);
            F_lig_j = f_stiffness(CSA,ls,l_lig);
    
            % calculate moment
            M_lig_j = F_lig_j.*d_lig;
            M_i = M_i + M_lig_j;
        end
    
        % sort angle and moment in order of angle
        [q_i,idx_sort] = sort(q_i);
        M_i = M_i(idx_sort);
    
        % create look-up table
        f_lig_single_coord = interpolant(['f_lig_single_coord_' model_info.ExtFunIO.coord_names.all{i}],...
            'bspline',{q_i},M_i);
        
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
        load nCoeffMat
        load expoVal_all
        for i=1:model_info.ligament_info.NLigament      
            index_dof_crossing = find(ligament_spanning_info_m(i,:)==1);
            nr_dof_crossing = length(index_dof_crossing);
            if nr_dof_crossing == 0
                continue
            end
            order = LigamentInfo_m.muscle(i).order;
            [mat,~] = n_art_mat_3_cas_SX_7(qin_SX(index_dof_crossing,1)',nCoeffMat(order,nr_dof_crossing),...
                expoVal_all{order,nr_dof_crossing});
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
                'the force-length of ligament ' model_info.ligament_info.parameters(j).name '.'])
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
    
    

else
    % if there are no ligaments, moments are zero (sparse)
    M_lig_multi_MX = MX(model_info.ExtFunIO.jointi.nq.all,1);
end

% total
M_lig_MX = M_lig_single_MX + M_lig_multi_MX;

f_ligamentMoment = Function('f_ligamentMoment',{Qs_MX},{M_lig_MX},{'Qs'},{'M_lig'});
f_ligamentMoment_single = Function('f_ligamentMoment_single',{Qs_MX},{M_lig_single_MX},{'Qs'},{'M_lig'});
f_ligamentMoment_multi = Function('f_ligamentMoment_multi',{Qs_MX},{M_lig_multi_MX},{'Qs'},{'M_lig'});

end % end of function