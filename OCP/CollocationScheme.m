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
%   * Collocation points
% 
%   - C -
%   * Coefficients of the collocation equation
% 
%   - D -
%   * Coefficients of the continuity equation
% 
%   - B -
%   * Coefficients of the quadrature function
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
% --------------------------------------------------------------------------
% This file is part of PredSim.
% 
% PredSim: A Framework for Rapid Predictive Simulations of Locomotion
% Copyright (c) 2026 KU Leuven
% 
% PredSim is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Affero General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% 
% PredSim is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public 
% License for more details.
% 
% You should have received a copy of the GNU Affero General Public License 
% along with PredSim. If not, see <https://www.gnu.org/licenses/>.
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
