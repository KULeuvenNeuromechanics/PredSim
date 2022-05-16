function [mat,diff_mat_q] = n_art_mat_3_cas_SX(q, order)
% --------------------------------------------------------------------------
% n_art_mat_3_cas_SX
%   This function returns the polynomials for approximating the muscle-tendon
%   lengths, velocities and moment-arms.
%   
% INPUT:
%   - q -
%   * 
%
%   - order -
%   * 
%
% OUTPUT:
%   - mat -
%   * 
%
%   - diff_mat_q -
%   * 
% 
% Original author: Original code from Wouter Aerts, adapted by Antoine Falisse
% Original date: 12/19/2018
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.*
n_dof = length(q(1,:));
q_all = SX(1,4);
for dof_nr = 1:n_dof
    q_all(1,dof_nr) = q(1,dof_nr);
end
temp_empty = 4-n_dof;
MX_null = SX(1);
MX_null(1,1) = 0;
for n_temp_empty = 1:temp_empty
    q_all(1,n_dof+n_temp_empty) = MX_null;
end

nr_coefficients = 0;
for n_q1 = 0:order
    if n_dof<2
        n_q2s = 0;
    else
        n_q2s = 0:order-n_q1;
    end
    for n_q2 = n_q2s
        if n_dof<3
            n_q3s = 0;
        else
            n_q3s = 0:order-n_q1-n_q2;
        end
        for n_q3 = n_q3s
            if n_dof<4
                n_q4s = 0;
            else
                n_q4s = 0:order-n_q1-n_q2-n_q3;
            end
            for n_q4 = n_q4s
                nr_coefficients = nr_coefficients + 1;
            end
        end
    end
end

coeff_nr = 1;
mat = SX(1,nr_coefficients);
diff_mat_q = SX(nr_coefficients,4);
for n_q1 = 0:order
    if n_dof<2
        n_q2s = 0;
    else
        n_q2s = 0:order-n_q1;
    end
    for n_q2 = n_q2s
        if n_dof<3
            n_q3s = 0;
        else
            n_q3s = 0:order-n_q1-n_q2;
        end
        for n_q3 = n_q3s
            if n_dof<4
                n_q4s = 0;
            else
                n_q4s = 0:order-n_q1-n_q2-n_q3;
            end           
            
            for n_q4 = n_q4s
                mat(1,coeff_nr) = q_all(:,1).^n_q1.*q_all(:,2).^n_q2.*q_all(:,3).^n_q3.*q_all(:,4).^n_q4;
                
                diff_mat_q1 = n_q1*q_all(:,1).^(n_q1-1).*q_all(:,2).^n_q2.*q_all(:,3).^n_q3.*q_all(:,4).^n_q4;
                diff_mat_q2 = q_all(:,1).^n_q1.*n_q2.*q_all(:,2).^(n_q2-1).*q_all(:,3).^n_q3.*q_all(:,4).^n_q4;
                diff_mat_q3 = q_all(:,1).^n_q1.*q_all(:,2).^n_q2.*n_q3.*q_all(:,3).^(n_q3-1).*q_all(:,4).^n_q4;
                diff_mat_q4 = q_all(:,1).^n_q1.*q_all(:,2).^n_q2.*q_all(:,3).^n_q3.*n_q4.*q_all(:,4).^(n_q4-1);
                
                diff_mat_q(coeff_nr,1) = diff_mat_q1;
                diff_mat_q(coeff_nr,2) = diff_mat_q2;
                diff_mat_q(coeff_nr,3) = diff_mat_q3;
                diff_mat_q(coeff_nr,4) = diff_mat_q4;                
               
                coeff_nr = coeff_nr + 1;
            end
        end
    end
end
end