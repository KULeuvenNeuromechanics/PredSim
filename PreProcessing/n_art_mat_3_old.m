% This function returns the polynomials for approximating the muscle-tendon
% lengths, velocities and moment-arms.
%
% Authors: Original code from Wouter Aerts, adapted by Antoine Falisse
% Date: 12/19/2018
%
function [mat,diff_mat_q] = n_art_mat_3_old(q, order)

n_dof = length(q(1,:));
nr_points = length(q(:,1));
q_all = zeros(nr_points, 4);
for dof_nr = 1:n_dof
    q_all(:,dof_nr) = q(:,dof_nr);
end
q = q_all;

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
mat = zeros(nr_points, nr_coefficients);
diff_mat_q = zeros(nr_points, nr_coefficients, n_dof);

coeff_nr = 1;
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
                mat(:,coeff_nr) = q(:,1).^n_q1.*q(:,2).^n_q2.*q(:,3).^n_q3.*q(:,4).^n_q4;
                
                diff_mat_q1 = n_q1*q(:,1).^(n_q1-1).*q(:,2).^n_q2.*q(:,3).^n_q3.*q(:,4).^n_q4;
                diff_mat_q2 = q(:,1).^n_q1.*n_q2.*q(:,2).^(n_q2-1).*q(:,3).^n_q3.*q(:,4).^n_q4;
                diff_mat_q3 = q(:,1).^n_q1.*q(:,2).^n_q2.*n_q3.*q(:,3).^(n_q3-1).*q(:,4).^n_q4;
                diff_mat_q4 = q(:,1).^n_q1.*q(:,2).^n_q2.*q(:,3).^n_q3.*n_q4.*q(:,4).^(n_q4-1);
                
                for dof_nr = 1:n_dof
                    eval(['diff_mat_q(:,coeff_nr,dof_nr) = diff_mat_q', num2str(dof_nr), ';']);
                end
                
                coeff_nr = coeff_nr + 1;
            end
        end
    end
end
end