function [struct_array] = double_array_to_struct_array(struct_array,fieldname,double_array)
% --------------------------------------------------------------------------
% double_array_to_struct_array
%   Helper function to read an array of doubles and add the values to a
%   field with given fieldname in the given cell array of structs.
% 
% INPUT:
%   - struct_array -
%   * Cell array of structs.
%
%   - fieldname -
%   * Name of a field in the structs. Should contain an array of doubles.
%
%   - double_array -
%   * Each row has the values from "fieldname" in the corresponding cell. 
%
% OUTPUT:
%   - struct_array -
%   * Cell array of structs. 
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


% number of structs in cell array
if isempty(struct_array)
    N = length(double_array);
else
    N = length(struct_array);
end

if isempty(double_array)
    for i=1:N
        struct_array(i).(fieldname) = [];
    end
else

    % check size mismatch
    if size(double_array,1)~=N
        if size(double_array,2)==N
            double_array = double_array';
        else
            error(['Size mismatch.'])
        end
    end

    for i=1:N
        array_i = double_array(i,:);
        if iscell(array_i) && numel(array_i)==1
            array_i = array_i{1};
        end
        struct_array(i).(fieldname) = array_i;
    end

end



