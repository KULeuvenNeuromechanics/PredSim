function [y,ydx] = mvpolyval(coeff, x, mu)
% --------------------------------------------------------------------------
% mvpolyval
%   Evaluate a multivariate polynomial.
%
%   The 2nd order polynomial in variables x1 and x2
%       y = c00 + c01*x2 + c02*x2^2 + c10*x1 + c11*x1*x2 + c20*x1^2
%   
%   can be written as the inner product of a vector of coefficients and a 
%   vector of monomials:
%       coeff = [c00; c01; c02; c10; c11; c20];
%       mno = [0; x2; x2^2; x1; x1*x2; x1^2];
%       y = mno' * coeff; (eq 1)
%
%   and its partial derivatives as
%       ydx = jac_mno' * coeff; (eq 2)
%
%   See also createMonimial mvpolyfit polyval
%
% INPUT:
%   - coeff -
%   * Vector with coefficients of the fitted polynimial.
%
%   - x -
%   * Independent variables. 
%     Concatenate vertically to evaluate multiple points.
%
%   - mu - (optional)
%   * Centering and scaling constants, obtained from mvpolyfit.
%
% OUTPUT:
%   - y -
%   * Polynomial evaluated at x
%
%   - ydx - (optional)
%   * Derivative of the polynomial, evaluated at x
% 
% Original author: Lars D'Hondt
% Original date: September 2025
% -------------------------------------------------------------------------
arguments
    coeff (:,1)
    x (:,:)
    mu (:,:) double = [0;1];
end

if (size(x,1)<length(coeff) && size(x,2)>=length(coeff))
    x = x';
end


x_sc = (x - mu(1,:)) ./mu(2,:);

order = find_order(length(coeff),size(x,2));

if nargout == 1
    mno = createMonomial(x_sc, order);

else % 1st derivative
    [mno,jac_mno] = createMonomial(x_sc, order, 2);
    % chain rule!
    jac_mno = jac_mno./repelem(mu(2,:)',numel(x)/size(mu,2),1);

    ydx = reshape(jac_mno*coeff,size(x));

end

y = mno*coeff;

end % end of function

function order = find_order(n_coeff, n_dof)
    p = 0;
    order = -1;
    lhs = factorial(n_dof) * n_coeff;
    while p < n_coeff
        rhs = prod(p+(1:n_dof));
        if abs(lhs - rhs) < 2*eps(lhs)
            order = p;
            break
        end
        p = p + 1;
    end
end
