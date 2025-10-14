function [mat,diff_mat_q] = n_art_mat(qs, order, diff_dim, diff_order)
% --------------------------------------------------------------------------
%
%
%   The vector of monomials can be expressed as a matrix of exponents
%       x1  x2
%        0   0
%        0   1
%        0   2
%        1   0
%        1   1
%        2   0
%
%
%
% --------------------------------------------------------------------------

if nargin == 2
    diff_dim = 3;
end
if nargin < 4
    diff_order = 1;
end
n_dof = size(qs,2);
n_points = size(qs,1);


monomial_exponents = generateCombinations(n_dof, order);

n_coeff = size(monomial_exponents,1);

if diff_dim == 0 || nargout == 1
    diff_mat_q = [];
end

if n_points == 1

    mat = eval_monomial(qs, monomial_exponents)';

    if diff_dim
        diff_mat_q = zeros(n_dof, n_coeff, 'like', qs);
        for j=1:n_dof
            diff_mat_q(j,:) = eval_der_monomial(qs, monomial_exponents, j, diff_order);
        end
    end

else
    mat = zeros(n_points, n_coeff, 'like', qs);
    for i=1:n_points
        mat(i,:) = eval_monomial(qs(i,:), monomial_exponents);
    end
    
    if diff_dim == 3
        diff_mat_q = zeros(n_points, n_coeff, n_dof, 'like', qs);
        for j=1:n_dof
            for i=1:n_points
                diff_mat_q(i,:,j) = eval_der_monomial(qs(i,:), monomial_exponents, j, diff_order);
            end
        end
    elseif diff_dim == 2
        diff_mat_q = zeros(n_points*n_dof, n_coeff, 'like', qs);
        for j=1:n_dof
            for i=1:n_points
                diff_mat_q(i+(j-1)*n_points,:) = eval_der_monomial(qs(i,:), monomial_exponents, j, diff_order);
            end
        end
    end
end

end

%%

% adapted from https://github.com/opensim-org/opensim-core/blob/main/
%   OpenSim/Common/MultivariatePolynomialFunction.cpp#L339
function combinations = generateCombinations(dimension, order)
    combinations = generateCombinationsHelper([], 1, 0, dimension, order, {});
    combinations = cell2mat(combinations');
end
function combinations = generateCombinationsHelper(current, level, sum,...
    dimension, order, combinations)
    if level > dimension
        if sum <= order
            combinations{end+1} = current;
        end
        return;
    end
    for i = 0:(order - sum)
        newCurrent = [current, i];
        combinations = generateCombinationsHelper(newCurrent, level + 1,...
            sum + i, dimension, order, combinations);
    end
end


function [y] = eval_monomial(x, powers)
    y = x(1) .^powers(:,1);
    for i=2:size(powers,2)
        y = y.* (x(i) .^powers(:,i));
    end
end

function [y] = eval_der_monomial(x, powers, idx, n)
    if nargin == 3
        n = 1;
    end
    coeff = ones(size(powers,1),1);
    for i=1:n
        coeff = coeff.*powers(:,idx);
        powers(:,idx) = powers(:,idx)-1;
        powers(powers<0) = 0;
    end
    y = x(1) .^powers(:,1).*coeff;
    for i=2:size(powers,2)
        y = y.* (x(i) .^powers(:,i));
    end
end