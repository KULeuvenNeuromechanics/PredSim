function [double_array] = struct_array_to_double_array(struct_array,fieldname)
% --------------------------------------------------------------------------
% struct_array_to_double_array
%   Helper function to read a cell array of structs and return an array
%   with the values of given fieldname.
% 
% INPUT:
%   - struct_array -
%   * Cell array of structs.
%
%   - fieldname -
%   * Name of a field in the structs. Should contain an array of doubles.
%
% OUTPUT:
%   - double_array -
%   * Each row has the values from "fieldname" in the corresponding cell. 
% 
% Original author: Lars D'Hondt
% Original date: 5/May/2022
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


N = length(struct_array);
Nval = 1; % expected number of doubles in field
double_array = nan(N,Nval);

for i=1:N
    if isfield(struct_array(i),fieldname) && ~isempty(struct_array(i).(fieldname))
        val_i = struct_array(i).(fieldname);
        Nval_i = numel(val_i);
        % add more columns if needed
        if Nval_i>Nval
            double_array = [double_array,nan(N,Nval_i-Nval)];
        end
        double_array(i,1:Nval_i) = reshape(val_i,[1,Nval_i]);
    end
end



