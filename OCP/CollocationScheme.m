function [tau_root,C,D,B] = CollocationScheme(d,method)
% --------------------------------------------------------------------------
% CollocationScheme
%    This script returns matrices needed for setting the collocation problem. 
%    This function is copied from an example (direct_collocation.m) from the 
%    CasADi example pack
%   
% INPUT:
%   - d -
%   * degree of the interpolating polynomial
%
%   - method -
%   * 'radau' or 'legendre'
%
% OUTPUT:
%   - tau_root -
%   * 
% 
%   - C -
%   * 
% 
%   - D -
%   * 
% 
%   - B -
%   * 
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.*

% Get collocation points
tau_root = [0 collocation_points(d,method)];

% Coefficients of the collocation equation
C = zeros(d+1,d+1);

% Coefficients of the continuity equation
D = zeros(d+1, 1);

% Coefficients of the quadrature function
B = zeros(d+1, 1);

% Construct polynomial basis
for j=1:d+1
    % Construct Lagrange polynomials to get the polynomial basis at the 
    % collocation point
    coeff = 1;
    for r=1:d+1
        if r ~= j
            coeff = conv(coeff, [1, -tau_root(r)]);
            coeff = coeff / (tau_root(j)-tau_root(r));
        end
    end
    % Evaluate the polynomial at the final time to get the coefficients of 
    % the continuity equation
    D(j) = polyval(coeff, 1.0);
    
    % Evaluate the time derivative of the polynomial at all collocation 
    % points to get the coefficients of the collocation equation
    pder = polyder(coeff);
    for r=1:d+1
        C(j,r) = polyval(pder, tau_root(r));
    end
    
    % Evaluate the integral of the polynomial to get the coefficients of 
    % the quadrature function
    pint = polyint(coeff);
    B(j) = polyval(pint, 1.0);
end
