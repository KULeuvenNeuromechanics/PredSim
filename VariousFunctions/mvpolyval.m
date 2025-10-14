function [y,ydx,varargout] = mvpolyval(coeff, x, mu)
arguments
    coeff (:,1) %double
    x (:,:)
    mu (:,:) double = [0;1];
end

if (size(x,1)<length(coeff) && size(x,2)>=length(coeff))
    x = x';
end


x_sc = (x - mu(1,:)) ./mu(2,:);

order = find_order(length(coeff),size(x,2));

if nargout >= 2
    [mno,jac_mno] = n_art_mat(x_sc, order, 2);
    % chain rule!
    jac_mno = jac_mno./repelem(mu(2,:)',numel(x)/size(mu,2),1);

    ydx = reshape(jac_mno*coeff,size(x));
else
    mno = n_art_mat(x_sc, order);
%     ydx = [];
end

if nargout >= 3
    [~,jac_mno] = n_art_mat(x_sc, order, 2, 3);
    % chain rule!
    jac_mno = jac_mno./repelem(mu(2,:)',numel(x)/size(mu,2),1);

    varargout{1} = reshape(jac_mno*coeff,size(x));
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

