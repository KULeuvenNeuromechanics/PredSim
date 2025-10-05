function [coeff, stats, mu] = mvpolyfit(X, Y, order, YdX, options)
% --------------------------------------------------------------------------
% mvpolyfit
%   Fit a multivariate polynomial.
%
%   The 2nd order polynomial in variables x1 and x2
%       y = c00 + c01*x2 + c02*x2^2 + c10*x1 + c11*x1*x2 + c20*x1^2
%   
%   Can be written as the inner product of a vector of coefficients and a 
%   vector of monomials:
%       coeff = [c00; c01; c02; c10; c11; c20];
%       mno = [0; x2; x2^2; x1; x1*x2; x1^2];
%       y = mno' * coeff;
%       
% 
%
%   References
%   [1] Harba and Serrancolí, 2025
%
%   See also n_art_mat mvpolyval polyfit
%
% INPUT:
%   - X -
%   * brief description of input_1, including data type.
% 
%   - Y -
%   * brief description of input_2. For more complex inputs, such as structs,
%   include an example input.
%
%   - order -
%   * 
%
%   - YdX - (optional) 
%   * brief description of input_3, including data type.
%
%   - options - (optional)
%   * name-value pairs with options
%
% OUTPUT:
%   - coeff -
%   * Vector with coefficients of the fitted polynimial.
%
%   - stats -
%   * Order of the fitted polynomial.
%
%   - mu -
%   * 
% 
% Original author: Lars D'Hondt
% Original date: September 2025
% --------------------------------------------------------------------------
arguments
    X (:,:) double
    Y (:,1) double
    order (1,:) {mustBeInteger, mustBeNonnegative}
    YdX (:,:) double = [];

    options.y_scale_mode {mustBeMember(options.y_scale_mode,...
        {'none','rms','max','mean','manual'})} = 'none';
    options.scale_y (1,1) double = 1;
    options.scale_ydx (1,:) double = ones(1,min(size(X)));

    options.threshold_rmse = 0.003;
    options.threshold_rmse_y = 0.003;
    options.threshold_rmse_ydx = 0.003;
    options.threshold_rel_err = 0.05;

    options.reduced_coeff = false;
    options.reduced_coeff_n_bins = 10;
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

%% scale data points
switch options.y_scale_mode
    case 'none'
        sf_y = 1;
        sf_ydx = ones(1,n_dof);
    case 'rms'
        sf_y = rms(Y);
        sf_ydx = rms(YdX,1);
    case 'max'
        sf_y = max(abs(Y));
        sf_ydx = max(abs(YdX),[],1);
    case 'mean'
        sf_y = abs(mean(Y));
        sf_ydx = abs(mean(YdX,1));
end


sf_jac = repelem(sf_ydx(:),n_points,1);

y_sc = Y./sf_y;
ydx_sc = YdX./sf_ydx;


if nargout >= 3
    mu = [mean(X,1); std(X,[],1)];
    x_sc = (X - mu(1,:)) ./mu(2,:);
else
    x_sc = X;
end

%% fit

fit_flag = 'max_order';

for fit_order=min(order):max(order)

    if isempty(YdX)
        [mno] = n_art_mat(x_sc, fit_order);
        X_aug = mno/sf_y;
        Y_aug = y_sc;
    else
        [mno,jac_mno] = n_art_mat(x_sc, fit_order, 2);
        if nargout >= 3
            % chain rule!
            jac_mno = jac_mno./repelem(mu(2,:)',n_points,1);
        end
        X_aug = [mno/sf_y; jac_mno./sf_jac];
        Y_aug = [y_sc; ydx_sc(:)];
    end
    
    coeff = X_aug \ Y_aug;
    
    
    rmse_aug = rms(Y_aug - X_aug*coeff);
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
        rel_ydx = [];
    end

    if (rmse_aug < options.threshold_rmse)
        fit_flag = 'rmse';
        break
    elseif ((rmse_y < options.threshold_rmse_y) && all(rmse_ydx < options.threshold_rmse_ydx))
        fit_flag = 'rmse_y';
        break
    elseif all(rel_y) && all(rel_ydx(:))
        fit_flag = 'rel_y';
        break
    end

end

%%

if nargout > 1
    stats.order = fit_order;
    stats.fit_flag = fit_flag;

    stats.sf_y = sf_y;
    stats.sf_ydx = sf_ydx;

    stats.rmse_aug = rmse_aug;
    stats.rmse_y = rmse_y;
    stats.rmse_ydx = rmse_ydx;
    stats.rel_err_y = err_y./abs(Y);
    stats.rel_err_ydx = err_ydx./abs(YdX);

    

end

%%

if options.reduced_coeff %&& strcmp(fit_flag,'rmse')
    % adapted from code by Mohanad Harba and Gil Serrancolí

    coeff_1 = coeff;

    threshold_rmse = max(rmse_aug,options.threshold_rmse);
    threshold_rmse_y = max(rmse_y,options.threshold_rmse_y);
    threshold_rmse_ydx = max(rmse_ydx,options.threshold_rmse_ydx);

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
        % add
        [bincount,~,coeff_idx] = histcounts(log10(p_values),...
            options.reduced_coeff_n_bins);
        idx_bins = find(bincount>0);
    end

    for i=[0,idx_bins]

        % Update selected coefficients
        selected_index = find(coeff_idx<=i);

        coeff_aux = X_aug(:, selected_index) \ Y_aug;
            % coefficients that best fit lMT and moment arms using only the 
            % subset of selected coefficients

        % Update the full coefficient vector with the selected ones
        coeff_2 = zeros(size(mno, 2), 1);
        coeff_2(selected_index) = coeff_aux;
        
        % Calculate rmse with the updated selected coefficients
        Y_reco_2 = X_aug * coeff_2;
        rmse_2 = rms(Y_aug-Y_reco_2);

        rmse_y_2 = rms(Y - mno*coeff_2);
        if ~isempty(YdX)
            rmse_ydx_2 = rms(YdX - reshape(jac_mno*coeff_2,size(YdX)),1);
        else
            rmse_ydx_2 = [];
        end

        if (rmse_2 < threshold_rmse) || ...
                ((rmse_y_2 < threshold_rmse_y) && all(rmse_ydx_2 < threshold_rmse_ydx))
            break
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
