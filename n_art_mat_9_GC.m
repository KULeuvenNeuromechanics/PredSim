function [mat,diff_mat_q] = n_art_mat_9_GC(q, order)

n_dof = length(q(1,:));
nr_points = length(q(:,1));
q_all = zeros(nr_points, 9);
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
                                    mat(:,coeff_nr) = q(:,1).^n_q1.*q(:,2).^n_q2.*q(:,3).^n_q3.*q(:,4).^n_q4.*q(:,5).^n_q5.*q(:,6).^n_q6.*q(:,7).^n_q7.*q(:,8).^n_q8.*q(:,9).^n_q9;

                                    diff_mat_q1 = n_q1*q(:,1).^(n_q1-1).*q(:,2).^n_q2.*q(:,3).^n_q3.*q(:,4).^n_q4.*q(:,5).^n_q5.*q(:,6).^n_q6.*q(:,7).^n_q7.*q(:,8).^n_q8.*q(:,9).^n_q9;
                                    diff_mat_q2 = q(:,1).^n_q1.*n_q2.*q(:,2).^(n_q2-1).*q(:,3).^n_q3.*q(:,4).^n_q4.*q(:,5).^n_q5.*q(:,6).^n_q6.*q(:,7).^n_q7.*q(:,8).^n_q8.*q(:,9).^n_q9;
                                    diff_mat_q3 = q(:,1).^n_q1.*q(:,2).^n_q2.*n_q3.*q(:,3).^(n_q3-1).*q(:,4).^n_q4.*q(:,5).^n_q5.*q(:,6).^n_q6.*q(:,7).^n_q7.*q(:,8).^n_q8.*q(:,9).^n_q9;
                                    diff_mat_q4 = q(:,1).^n_q1.*q(:,2).^n_q2.*q(:,3).^n_q3.*n_q4.*q(:,4).^(n_q4-1).*q(:,5).^n_q5.*q(:,6).^n_q6.*q(:,7).^n_q7.*q(:,8).^n_q8.*q(:,9).^n_q9;
                                    diff_mat_q5 = q(:,1).^n_q1.*q(:,2).^n_q2.*q(:,3).^n_q3.*q(:,4).^n_q4.*n_q5.*q(:,5).^(n_q5-1).*q(:,6).^n_q6.*q(:,7).^n_q7.*q(:,8).^n_q8.*q(:,9).^n_q9;
                                    diff_mat_q6 = q(:,1).^n_q1.*q(:,2).^n_q2.*q(:,3).^n_q3.*q(:,4).^n_q4.*q(:,5).^n_q5.*n_q6.*q(:,6).^(n_q6-1).*q(:,7).^n_q7.*q(:,8).^n_q8.*q(:,9).^n_q9;
                                    diff_mat_q7 = q(:,1).^n_q1.*q(:,2).^n_q2.*q(:,3).^n_q3.*q(:,4).^n_q4.*q(:,5).^n_q5.*q(:,6).^n_q6.*n_q7.*q(:,7).^(n_q7-1).*q(:,8).^n_q8.*q(:,9).^n_q9;
                                    diff_mat_q8 = q(:,1).^n_q1.*q(:,2).^n_q2.*q(:,3).^n_q3.*q(:,4).^n_q4.*q(:,5).^n_q5.*q(:,6).^n_q6.*q(:,7).^n_q7.*n_q8.*q(:,8).^(n_q8-1).*q(:,9).^n_q9;
                                    diff_mat_q9 = q(:,1).^n_q1.*q(:,2).^n_q2.*q(:,3).^n_q3.*q(:,4).^n_q4.*q(:,5).^n_q5.*q(:,6).^n_q6.*q(:,7).^n_q7.*q(:,8).^n_q8.*n_q9.*q(:,9).^(n_q9-1);

                                    for dof_nr = 1:n_dof
                                        eval(['diff_mat_q(:,coeff_nr,dof_nr) = diff_mat_q', num2str(dof_nr),';']);
                                    end
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
