function [nr_coefficients,expoVal]=calcCoeff(order,n_dof)
% --------------------------------------------------------------------------
% calcCoeff 
% Function to compute the helper matricies used in making CasADi functio
% that determine the muscle tendon length, velocity and moment arm using
% polynomial fit.
% 
% INPUT:
%   - order -
%   * order of fit
% 
%   - n_dof -
%   * number of degrees of freedom
% 
% OUTPUT:
%   - nr_coefficients -
%   * number of coefficients needed.
% 
%   - expoVal -
%   * valuse of exponents.
% 
% Original authors: Dhruv Gupta
% Original date: 30/11/2021
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

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
                if n_dof<5
                    n_q5s = 0;
                else
                    n_q5s = 0:order-n_q1-n_q2-n_q3-n_q4;
                end
                for n_q5 = n_q5s
                    if n_dof<6
                        n_q6s = 0;
                    else
                        n_q6s = 0:order-n_q1-n_q2-n_q3-n_q4-n_q5;
                    end
                    for n_q6 = n_q6s
                        nr_coefficients = nr_coefficients + 1;
                        expoVal(nr_coefficients,1) = n_q1;
                        expoVal(nr_coefficients,2) = n_q2;
                        expoVal(nr_coefficients,3) = n_q3;
                        expoVal(nr_coefficients,4) = n_q4;
                        expoVal(nr_coefficients,5) = n_q5;
                        expoVal(nr_coefficients,6) = n_q6;
                    end
                end
            end
        end
    end
end
