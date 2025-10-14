function [coeff, stats, mu] = mvpolyfit(X, Y, order, YdX, options)
% --------------------------------------------------------------------------
% mvpolyfit
%   Fit a multivariate polynomial.
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
%   Given data samples X, Y, coeff can be estimated from a least-squares
%   fit w.r.t. eq 1. If derivative information is known, eq 2 can be
%   included to improve the fit.
%
%   Terms that do not contribute to the goodness of fit can be excluded [1].
% 
%
%   References
%   [1] Harba and Serrancolí, 2025
%
%   See also n_art_mat mvpolyval polyfit
%
% INPUT:
%   - X -
%   * Matrix with sampled values of independent variables.
% 
%   - Y -
%   * Vector with sampled values of dependent variable.
%
%   - order -
%   * Order of the polynomial to be fitted. 
%
%   - YdX - (optional) 
%   * Matrix with sampled values of the Jacobian of Y to X.
%
%   - options - (optional)
%   * name-value pairs with options
%       - threshold_rmse_y      rmse on eq 1
%       - threshold_rmse_ydx    rmse on eq 2
%       - threshold_rel_err     relative error on each samples
%       - reduced_coeff         try to reduce the amount of coefficients
%       - reduced_coeff_n_bins  divide the coefficients in bins instead of
%                               adding one by one when trying to reduce the 
%                               amount of coefficients.
%       - mldivide              custom function to perform A\b
%
% OUTPUT:
%   - coeff -
%   * Vector with coefficients of the fitted polynimial.
%
%   - stats - (optional)
%   * Struct with information about the fit.
%
%   - mu - (optional)
%   * Perform centering and scaling of X before fitting, and return the
%   centering and scaling values. If this is used, mu needs to be passed to
%   mvpolyval to correctly evaluate the fitted polynomial.
%   Centering and scaling are chosen such that the resulting X has
%   zero-mean and unit standard deviation, as this improves the numerical
%   stability.
% 
% Original author: Lars D'Hondt
% Original date: September 2025
% --------------------------------------------------------------------------
arguments
    X (:,:) double
    Y (:,1) double
    order (1,:) {mustBeInteger, mustBeNonnegative}
    YdX (:,:) double = [];

    % if order is a range of orders
    options.threshold_rmse_y = 0.003;
    options.threshold_rmse_ydx = 0.003;
    options.threshold_rel_err = 0.05;

    options.reduced_coeff = false;
    options.reduced_coeff_n_bins = 10;

    options.mldivide = [];
end

%% check input sizes
n_points = length(Y);

if size(X,2) == n_points
    X = X';
elseif size(X,1) ~= n_points
    error("size mismatch X and y")
end

n_dof = size(X,2);

if ~isempty(YdX)
    if size(YdX,1) == n_dof && size(YdX,2) == n_points
        YdX = YdX';
    elseif ~(size(YdX,1) == n_points && size(YdX,2) == n_dof)
        error("size mismatch ydx with x and y")
    end
end

if ~isempty(options.mldivide)
    if ischar(options.mldivide) || isstring(options.mldivide)
        mldivide_impl = str2func(options.mldivide);
    elseif isa(options.mldivide,'function_handle')
        mldivide_impl = options.mldivide;
    end
end

%% scale data points

if nargout >= 3
    mu = [mean(X,1); std(X,[],1)];
    x_sc = (X - mu(1,:)) ./mu(2,:);
else
    x_sc = X;
end


%% fit

fit_flag = 'max_order';

for fit_order=min(order):max(order)

    % Create matrices
    if isempty(YdX)
        [mno] = n_art_mat(x_sc, fit_order);
        X_aug = mno;
        Y_aug = Y;
    else
        [mno,jac_mno] = n_art_mat(x_sc, fit_order, 2);
        if nargout >= 3
            % chain rule!
            jac_mno = jac_mno./repelem(mu(2,:)',n_points,1);
        end
        X_aug = [mno; jac_mno];
        Y_aug = [Y; YdX(:)];
    end
    
    % Fit coeff
    if isempty(options.mldivide)
        coeff = X_aug \ Y_aug;
    else
        coeff = mldivide_impl(X_aug, Y_aug);
    end
    
    % Evaluate acceptance criteria
    Y_fit = mno*coeff;
    rmse_y = rms(Y - Y_fit);
    err_y = abs(Y - Y_fit);
    rel_y = err_y < options.threshold_rel_err*abs(Y);
    if ~isempty(YdX)
        YdX_fit = reshape(jac_mno*coeff,size(YdX));
        rmse_ydx = rms(YdX - YdX_fit,1);
        err_ydx = abs(YdX - YdX_fit);
        rel_ydx = err_ydx < options.threshold_rel_err*abs(YdX);
    else
        rmse_ydx = [];
        err_ydx = [];
        rel_ydx = [];
    end

    if ((rmse_y < options.threshold_rmse_y) && all(rmse_ydx < options.threshold_rmse_ydx))
        fit_flag = 'rmse';
        break
    elseif all(rel_y) && all(rel_ydx(:))
        fit_flag = 'rel_err';
        break
    end

end

%%

if nargout > 1
    stats.order = fit_order;
    stats.fit_flag = fit_flag;

    stats.rmse_y = rmse_y;
    stats.rmse_ydx = rmse_ydx;
    stats.rel_err_y = err_y./abs(Y);
    stats.rel_err_ydx = err_ydx./abs(YdX);

end

%%

if options.reduced_coeff
    % adapted from code by Mohanad Harba and Gil Serrancolí

    coeff_1 = coeff;

    % Compute Y_hat
    y_hat = X_aug * coeff;
    SSE = sum((Y_aug - y_hat).^2);
    
    % Compute standard error (SE)
    jacobian_1 = pinv(X_aug' * X_aug);
    SE = sqrt(diag(SSE / (length(Y_aug) - length(coeff)) * jacobian_1)); 

    % Calculate the t-values for each coefficient
    t_values = coeff ./ SE;
    
    % Calculate the p-values for each coefficient using the t-distribution
    df = length(X_aug) - length(coeff);
    p_values = 2 * (1 - tcdf(abs(t_values), df));
    

    % Add coefficients (in the order of p-value) until the threshold is
    % satisfied.
    if options.reduced_coeff_n_bins <= 0
        % add one by one
        [~, index] = sort(p_values);
        coeff_idx = index-1;
        idx_bins = 1:length(coeff)-1;
    else
        % add multiple
        [bincount,~,coeff_idx] = histcounts(log10(p_values),...
            options.reduced_coeff_n_bins);
        idx_bins = find(bincount>0);
    end

    for i=[0,idx_bins]

        % Update selected coefficients
        selected_index = find(coeff_idx<=i);

        % Coefficient values that best fit data using only the subset of 
        % selected coefficients
        if isempty(options.mldivide)
            coeff_aux = X_aug(:, selected_index) \ Y_aug;
        else
            coeff_aux = mldivide_impl(X_aug(:, selected_index), Y_aug);
        end
            

        % Update the full coefficient vector with the selected ones
        coeff_2 = zeros(size(mno, 2), 1);
        coeff_2(selected_index) = coeff_aux;
        
        % Check acceptance criterium with updated selected coefficients.
        Y_fit = mno*coeff_2;

        switch fit_flag
            case 'rmse_y'
                rmse_y = rms(Y - Y_fit);
                if ~isempty(YdX)
                    YdX_fit = reshape(jac_mno*coeff_2,size(YdX));
                    rmse_ydx = rms(YdX - YdX_fit,1);
                else
                    rmse_ydx = [];
                end
                if ((rmse_y < options.threshold_rmse_y) && ...
                        all(rmse_ydx < options.threshold_rmse_ydx))
                    break
                end

            case 'rel_y'
                err_y = abs(Y - Y_fit);
                rel_y = err_y < options.threshold_rel_err*abs(Y);
                if ~isempty(YdX)
                    YdX_fit = reshape(jac_mno*coeff_2,size(YdX));
                    err_ydx = abs(YdX - YdX_fit);
                    rel_ydx = err_ydx < options.threshold_rel_err*abs(YdX);
                else
                    rel_ydx = [];
                end
                if all(rel_y) && all(rel_ydx(:))
                    break
                end

        end
        
    end

    if nargout > 1
        stats.coeff_full = coeff_1;
        stats.coeff_p_values = p_values;
        stats.coeff_p_bin = coeff_idx;
        stats.coeff_used = sort(selected_index);

        stats.rmse_red_aug = rmse_2;
        stats.rmse_red_y = rmse_y_2;
        stats.rmse_red_ydx = rmse_ydx_2;
    end
    coeff = coeff_2;

end

end % end of function
