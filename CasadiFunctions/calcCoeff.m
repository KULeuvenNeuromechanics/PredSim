function [nr_coefficients,expoVal]=calcCoeff(order,n_dof)
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
                expoVal(nr_coefficients,1) = n_q1;
                expoVal(nr_coefficients,2) = n_q2;
                expoVal(nr_coefficients,3) = n_q3;
                expoVal(nr_coefficients,4) = n_q4;
            end
        end
    end
end
