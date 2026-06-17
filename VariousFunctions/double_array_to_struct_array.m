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
% Last edit by: 
% Last edit date: 
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
        struct_array(i).(fieldname) = double_array(i,:);
    end

end



