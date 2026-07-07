function [F] = load_external_function(S)
% --------------------------------------------------------------------------
% load_external_function
%   The external function performs inverse dynamics through the
%   OpenSim/Simbody C++ API. This external function is compiled as a dll from
%   which we create a Function instance using CasADi in MATLAB.
% 
% INPUT:
%   - S -
%   * setting structure S
%
% OUTPUT:
%   - F -
%   * function handle to the external function
% 
% Original author: Lars D'Hondt
% Original date: 10/May/2022
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
pathmain = pwd;
cd(S.misc.subject_path)
F  = external('F',S.misc.external_function);
cd(pathmain);

