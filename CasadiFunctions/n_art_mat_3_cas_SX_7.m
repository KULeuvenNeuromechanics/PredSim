function [mat,diff_mat_q] = n_art_mat_3_cas_SX_7(q, nr_coefficients,expoVal)
% This function returns the polynomials for approximating the muscle-tendon
% lengths, velocities and moment-arms.
%
% Authors: Original code from Wouter Aerts, adapted by Antoine Falisse
% Date: 12/19/2018
% Simplified and changed to be able to handle 6 dof by Dhruv Gupta (May 12, 2022)

import casadi.*
n_dof = length(q(1,:));
q_all = SX(1,6);
q_all(1,1:n_dof) = q;
% if n_dof<4
%     temp_empty = 4-n_dof;
%     MX_null = MX(1,temp_empty);
%     MX_null(1,1:temp_empty) = 0;
%     q_all=[q_all,MX_null];
% end
temp_empty = 6-n_dof;
MX_null = SX(1);
MX_null(1,1) = 0;
for n_temp_empty = 1:temp_empty
    q_all(1,n_dof+n_temp_empty) = MX_null;
end

mat = SX(nr_coefficients,1);
diff_mat_q = SX(nr_coefficients,6);

q_all_nr = repmat(q_all,nr_coefficients,1);
expo_q_all_nr = q_all_nr.^expoVal;
mat(1:nr_coefficients,1) = expo_q_all_nr(:,1).*expo_q_all_nr(:,2).*expo_q_all_nr(:,3).*expo_q_all_nr(:,4).*expo_q_all_nr(:,5).*expo_q_all_nr(:,6);
diff_mat_q(:,1)=expo_q_all_nr(:,2).*expo_q_all_nr(:,3).*expo_q_all_nr(:,4).*expo_q_all_nr(:,5).*expo_q_all_nr(:,6).*expoVal(:,1).*(q_all_nr(:,1).^(expoVal(:,1)-1));
diff_mat_q(:,2)=expo_q_all_nr(:,1).*expo_q_all_nr(:,3).*expo_q_all_nr(:,4).*expo_q_all_nr(:,5).*expo_q_all_nr(:,6).*expoVal(:,2).*(q_all_nr(:,2).^(expoVal(:,2)-1));
diff_mat_q(:,3)=expo_q_all_nr(:,1).*expo_q_all_nr(:,2).*expo_q_all_nr(:,4).*expo_q_all_nr(:,5).*expo_q_all_nr(:,6).*expoVal(:,3).*(q_all_nr(:,3).^(expoVal(:,3)-1));
diff_mat_q(:,4)=expo_q_all_nr(:,1).*expo_q_all_nr(:,2).*expo_q_all_nr(:,3).*expo_q_all_nr(:,5).*expo_q_all_nr(:,6).*expoVal(:,4).*(q_all_nr(:,4).^(expoVal(:,4)-1));
diff_mat_q(:,5)=expo_q_all_nr(:,1).*expo_q_all_nr(:,2).*expo_q_all_nr(:,3).*expo_q_all_nr(:,4).*expo_q_all_nr(:,6).*expoVal(:,5).*(q_all_nr(:,5).^(expoVal(:,5)-1));
diff_mat_q(:,6)=expo_q_all_nr(:,1).*expo_q_all_nr(:,2).*expo_q_all_nr(:,3).*expo_q_all_nr(:,4).*expo_q_all_nr(:,5).*expoVal(:,6).*(q_all_nr(:,6).^(expoVal(:,6)-1));

end
