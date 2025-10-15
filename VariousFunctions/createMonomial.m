function [mno, der_mno] = createMonomial(x, order, der_dim, der_order)
% --------------------------------------------------------------------------
% createMonomial
%   Create a vector (matrix) of monomials
%
%   The 2nd order polynomial in variables x1 and x2
%       y = c00 + c01*x2 + c02*x2^2 + c10*x1 + c11*x1*x2 + c20*x1^2
%   
%   can be written as the inner product of a vector of coefficients and a 
%   vector of monomials:
%       coeff = [c00; c01; c02; c10; c11; c20];
%       mno = [0; x2; x2^2; x1; x1*x2; x1^2];
%       y = mno' * coeff;
%       
%   The vector of monomials can be expressed as a matrix of powers
%       x1  x2
%        0   0
%        0   1
%        0   2
%        1   0
%        1   1
%        2   0
%
%
%   Derivatives can be calculated from the chain rule:
%       ydx = jac_mno' * coeff;
%   With
%       jac_mno = [dmno_dx1, dmno_dx2]
%               = [0,    0; 
%                  0,    1; 
%                  0,  2*x2; 
%                  1,    0; 
%                  x2,  x1; 
%                  2*x1, 0];
% 
%
%   See also mvpolyfit mvpolyval
%
% INPUT:
%   - x -
%   * Horizontal vector of independent variables. 
%     Concatenate vertically to evaluate multiple.
% 
%   - order -
%   * Maximum order of the monomials. 
%
%   - der_dim - (optional) Default: 2
%   * Dimensions of the matrix with derivatives. 
%       2: Derivatives matrices for samples are concatenated vertically.
%       3: Derivatives matrices for samples are concatenated in the 3rd dim.
%
%   - der_order - (optional) Default: 1
%   * Matrix with sampled values of the Jacobian of Y to X.
%
% OUTPUT:
%   - mno -
%   * Vector (matrix) of monomials.
%
%   - der_mno - (optional)
%   * Struct with information about the fit.
% 
% Original author: Lars D'Hondt
% Original date: September 2025
% --------------------------------------------------------------------------
arguments
    x (:,:)
    order (1,1) {mustBeInteger, mustBeNonnegative}
    der_dim (1,1) {mustBeInteger, mustBeNonnegative} = 2;
    der_order (1,1) {mustBeInteger, mustBeNonnegative} = 1;
end

if size(x,2)==1 && size(x,1)>=1
    x = x';
end

n_dof = size(x,2);
n_points = size(x,1);




%% Matrix with powers
% Generate all combinations of 0:order with sum<=order
% adapted from https://github.com/opensim-org/opensim-core/blob/main/
%   OpenSim/Common/MultivariatePolynomialFunction.cpp
mno_pow = generate_combinations([], 1, 0, n_dof, order, {});
mno_pow = cell2mat(mno_pow');

n_coeff = size(mno_pow,1);

%% Monomial
mno = zeros(n_points, n_coeff, 'like', x);
for i=1:n_points
    mno(i,:) = eval_monomial(x(i,:), mno_pow);
end

%% Derivative
if nargout == 1
    return
end

if der_dim == 2
    der_mno = zeros(n_points*n_dof, n_coeff, 'like', x);
    for j=1:n_dof
        for i=1:n_points
            der_mno(i+(j-1)*n_points,:) = ...
                eval_der_monomial(x(i,:), mno_pow, j, der_order);
        end
    end

elseif der_dim == 3
    der_mno = zeros(n_points, n_coeff, n_dof, 'like', x);
    for j=1:n_dof
        for i=1:n_points
            der_mno(i,:,j) = ...
                eval_der_monomial(x(i,:), mno_pow, j, der_order);
        end
    end

end

end % end of function


%% helper functions

function combinations = generate_combinations(current, level, sum,...
    dimension, order, combinations)
    if level > dimension
        if sum <= order
            combinations{end+1} = current;
        end
        return;
    end
    for i = 0:(order - sum)
        newCurrent = [current, i];
        combinations = generate_combinations(newCurrent, level + 1,...
            sum + i, dimension, order, combinations);
    end
end


function [y] = eval_monomial(x, powers)
    y = x(1) .^powers(:,1);
    for i=2:size(powers,2)
        y = y.* (x(i) .^powers(:,i));
    end
end

function [y] = eval_der_monomial(x, powers, idx, der)
    if nargin == 3
        der = 1;
    end
    coeff = ones(size(powers,1),1);
    for i=1:der
        % y = c*x^p
        % dydx = p*c*x^(p-1)
        coeff = coeff.*powers(:,idx);
        powers(:,idx) = powers(:,idx)-1;
        powers(powers<0) = 0;
    end

    y = eval_monomial(x,powers).*coeff;
end
