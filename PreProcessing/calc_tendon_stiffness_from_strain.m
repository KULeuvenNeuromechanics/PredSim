function [kT] = calc_tendon_stiffness_from_strain(tendon_strain_at_FMo)


import casadi.*


lTtilde = tendon_strain_at_FMo + 1;

kT = SX.sym('kT',1,1);
shift = getShift(kT);

fse = (exp(kT.*(lTtilde - 0.995)))/5 - 0.25 + shift;

f_fse = Function('f_fse',{kT},{fse-1});

f_kT = rootfinder('f_kT','newton',f_fse);

kT = full(f_kT(35));

% if round(kT) == 35
%     kT = 35;
% end

end % end of function