function [mat,diff_mat_q] = n_art_mat_9_GC_SX(q, order)
import casadi.*
n_dof = length(q(1,:));
q_all = SX(1, 9);
for dof_nr = 1:n_dof
    q_all(1,dof_nr) = q(1,dof_nr);
end
temp_empty = 9-n_dof;
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
                        if n_dof<7
                            n_q7s = 0;
                        else
                            n_q7s = 0:order-n_q1-n_q2-n_q3-n_q4-n_q5-n_q6;
                        end
                        for n_q7 = n_q7s
                            if n_dof<8
                                n_q8s = 0;
                            else
                                n_q8s = 0:order-n_q1-n_q2-n_q3-n_q4-n_q5-n_q6-n_q7;
                            end
                            for n_q8 = n_q8s
                                if n_dof<9
                                    n_q9s = 0;
                                else
                                    n_q9s = 0:order-n_q1-n_q2-n_q3-n_q4-n_q5-n_q6-n_q7-n_q8;
                                end
                                for n_q9 = n_q9s
                                    nr_coefficients = nr_coefficients + 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% This is how Gil had it implemetned.
coeff_nr = 1;
mat = SX(1, nr_coefficients);
diff_mat_q = SX(nr_coefficients, 9);

% changed to try adn amtch old
%nr_points = 1;
%mat = zeros(nr_points, nr_coefficients);
%diff_mat_q = zeros(nr_points, nr_coefficients, n_dof);
%coeff_nr = 1;

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
                        if n_dof<7
                            n_q7s = 0;
                        else
                            n_q7s = 0:order-n_q1-n_q2-n_q3-n_q4-n_q5-n_q6;
                        end
                        for n_q7 = n_q7s
                            if n_dof<8
                                n_q8s = 0;
                            else
                                n_q8s = 0:order-n_q1-n_q2-n_q3-n_q4-n_q5-n_q6-n_q7;
                            end
                            for n_q8 = n_q8s
                                if n_dof<9
                                    n_q9s = 0;
                                else
                                    n_q9s = 0:order-n_q1-n_q2-n_q3-n_q4-n_q5-n_q6-n_q7-n_q8;
                                end
                                for n_q9 = n_q9s
                                    mat(1,coeff_nr) = q_all(:,1).^n_q1.*q_all(:,2).^n_q2.*q_all(:,3).^n_q3.*q_all(:,4).^n_q4.*q_all(:,5).^n_q5.*q_all(:,6).^n_q6.*q_all(:,7).^n_q7.*q_all(:,8).^n_q8.*q_all(:,9).^n_q9;

                                    diff_mat_q1 = n_q1*q_all(:,1).^(n_q1-1).*q_all(:,2).^n_q2.*q_all(:,3).^n_q3.*q_all(:,4).^n_q4.*q_all(:,5).^n_q5.*q_all(:,6).^n_q6.*q_all(:,7).^n_q7.*q_all(:,8).^n_q8.*q_all(:,9).^n_q9;
                                    diff_mat_q2 = q_all(:,1).^n_q1.*n_q2.*q_all(:,2).^(n_q2-1).*q_all(:,3).^n_q3.*q_all(:,4).^n_q4.*q_all(:,5).^n_q5.*q_all(:,6).^n_q6.*q_all(:,7).^n_q7.*q_all(:,8).^n_q8.*q_all(:,9).^n_q9;
                                    diff_mat_q3 = q_all(:,1).^n_q1.*q_all(:,2).^n_q2.*n_q3.*q_all(:,3).^(n_q3-1).*q_all(:,4).^n_q4.*q_all(:,5).^n_q5.*q_all(:,6).^n_q6.*q_all(:,7).^n_q7.*q_all(:,8).^n_q8.*q_all(:,9).^n_q9;
                                    diff_mat_q4 = q_all(:,1).^n_q1.*q_all(:,2).^n_q2.*q_all(:,3).^n_q3.*n_q4.*q_all(:,4).^(n_q4-1).*q_all(:,5).^n_q5.*q_all(:,6).^n_q6.*q_all(:,7).^n_q7.*q_all(:,8).^n_q8.*q_all(:,9).^n_q9;
                                    diff_mat_q5 = q_all(:,1).^n_q1.*q_all(:,2).^n_q2.*q_all(:,3).^n_q3.*q_all(:,4).^n_q4.*n_q5.*q_all(:,5).^(n_q5-1).*q_all(:,6).^n_q6.*q_all(:,7).^n_q7.*q_all(:,8).^n_q8.*q_all(:,9).^n_q9;
                                    diff_mat_q6 = q_all(:,1).^n_q1.*q_all(:,2).^n_q2.*q_all(:,3).^n_q3.*q_all(:,4).^n_q4.*q_all(:,5).^n_q5.*n_q6.*q_all(:,6).^(n_q6-1).*q_all(:,7).^n_q7.*q_all(:,8).^n_q8.*q_all(:,9).^n_q9;
                                    diff_mat_q7 = q_all(:,1).^n_q1.*q_all(:,2).^n_q2.*q_all(:,3).^n_q3.*q_all(:,4).^n_q4.*q_all(:,5).^n_q5.*q_all(:,6).^n_q6.*n_q7.*q_all(:,7).^(n_q7-1).*q_all(:,8).^n_q8.*q_all(:,9).^n_q9;
                                    diff_mat_q8 = q_all(:,1).^n_q1.*q_all(:,2).^n_q2.*q_all(:,3).^n_q3.*q_all(:,4).^n_q4.*q_all(:,5).^n_q5.*q_all(:,6).^n_q6.*q_all(:,7).^n_q7.*n_q8.*q_all(:,8).^(n_q8-1).*q_all(:,9).^n_q9;
                                    diff_mat_q9 = q_all(:,1).^n_q1.*q_all(:,2).^n_q2.*q_all(:,3).^n_q3.*q_all(:,4).^n_q4.*q_all(:,5).^n_q5.*q_all(:,6).^n_q6.*q_all(:,7).^n_q7.*q_all(:,8).^n_q8.*n_q9.*q_all(:,9).^(n_q9-1);

                                    diff_mat_q(coeff_nr,1) = diff_mat_q1;
                                    diff_mat_q(coeff_nr,2) = diff_mat_q2;
                                    diff_mat_q(coeff_nr,3) = diff_mat_q3;
                                    diff_mat_q(coeff_nr,4) = diff_mat_q4; 
                                    diff_mat_q(coeff_nr,5) = diff_mat_q5;
                                    diff_mat_q(coeff_nr,6) = diff_mat_q6;
                                    diff_mat_q(coeff_nr,7) = diff_mat_q7;
                                    diff_mat_q(coeff_nr,8) = diff_mat_q8;
                                    diff_mat_q(coeff_nr,9) = diff_mat_q9;
                                    
                                    coeff_nr = coeff_nr + 1;

                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
end


% function [mat, diff_mat_q] = n_art_mat_9_GC(q, order)
% 
% n_dof = length(q(1,:));
% nr_points = length(q(:,1));
% 
% % Calculate the number of coefficients
% nr_coefficients = nchoosek(n_dof + order, order);
% 
% % Initialize the mat and diff_mat_q matrices
% mat = zeros(nr_points, nr_coefficients);
% diff_mat_q = zeros(nr_points, nr_coefficients, n_dof);
% coeff_nr = 1;
% 
% % Loop over the degrees of freedom
% for n_q1 = 0:order
%     for n_q2 = 0:order - n_q1
%         for n_q3 = 0:order - n_q1 - n_q2
%             for n_q4 = 0:order - n_q1 - n_q2 - n_q3
%                 for n_q5 = 0:order - n_q1 - n_q2 - n_q3 - n_q4
%                     for n_q6 = 0:order - n_q1 - n_q2 - n_q3 - n_q4 - n_q5
%                         for n_q7 = 0:order - n_q1 - n_q2 - n_q3 - n_q4 - n_q5 - n_q6
%                             for n_q8 = 0:order - n_q1 - n_q2 - n_q3 - n_q4 - n_q5 - n_q6 - n_q7
%                                 for n_q9 = 0:order - n_q1 - n_q2 - n_q3 - n_q4 - n_q5 - n_q6 - n_q7 - n_q8
%                                     mat(:, coeff_nr) = q(:, 1).^n_q1 .* q(:, 2).^n_q2 .* q(:, 3).^n_q3 .* ...
%                                     q(:, 4).^n_q4 .* q(:, 5).^n_q5 .* q(:, 6).^n_q6 .* ...
%                                     q(:, 7).^n_q7 .* q(:, 8).^n_q8 .* q(:, 9).^n_q9;
% 
%                                     % Calculate the partial derivatives for each degree of freedom
%                                     for dof_nr = 1:n_dof
%                                         eval(['diff_mat_q(:, coeff_nr, dof_nr) = n_q1 * q(:, 1).^(n_q1 - 1) .* q(:, 2).^n_q2 .*'...
%                                               'q(:, 3).^n_q3 .* q(:, 4).^n_q4 .* q(:, 5).^n_q5 .* q(:, 6).^n_q6 .*'...
%                                               'q(:, 7).^n_q7 .* q(:, 8).^n_q8 .* q(:, 9).^n_q9;']);
%                                     end
% 
%                                     coeff_nr = coeff_nr + 1;
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% end
