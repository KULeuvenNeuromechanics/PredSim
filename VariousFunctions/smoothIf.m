function y = smoothIf(x,threshold_true,threshold_false)
% --------------------------------------------------------------------------
% smoothIf
%   Smooth approximation of if-statement, that is compatible with
%   algorithmic differentiation.
%
%   The original non-smooth statement
%       if x >= threshold
%           y = 1; % true
%       else
%           y = 0; % false
%       end
%
%   is written as
%       if x >= threshold_true
%           y = 1; % true
%       elseif x <= threshold_false
%           y = 0; % false
%       else
%           y > 0 & y < 1; % transient
%       end
%
%   To get true if x <= threshold, set threshold_true < threshold_false.
%   Both thresholds cannot be equal.
%
% INPUT:
%   - x -
%   * value to be compared to threshold
%
%   - threshold_true -
%   * threshold to consider the statement true
%
%   - threshold_false -
%   * threshold to consider the statement false
%
% OUTPUT:
%   - y -
%   * equals 1 if the statement is true, 0 if it is false
%
% 
% Original author: Lars D'Hondt
% Original date: 4/April/2023
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

% range that is affected by smoothing
range = threshold_true - threshold_false;

% midpoint of  range
mid = (threshold_false + threshold_true)/2;

% scale x to range
x_scaled = (x - mid)/range + 1/2;

% apply tanh-smoothing
y = (tanh((2*x_scaled-1)*pi)+1)/2;

end