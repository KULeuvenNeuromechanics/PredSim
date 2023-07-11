function [varargout] = unpack_name_value_combinations(name_values_cell, valid_names, value_sizes)
% --------------------------------------------------------------------------
% unpack_name_value_combinations
%   Helper function to unpack a cell array of name-value pairs,
%   name-value-value triplets, etc. Each set of values is returned as an
%   array with a row for each valid name. Unspecified values are NaN.
%   Unit tests in \Tests\test_unpack_name_value_combinations.m show how to
%   call this function.
% 
% INPUT:
%   - name_values_cell -
%   * cell array of name-value pairs, name-value-value triplets, etc.
%
%   - valid_names -
%   * cell array of names that are accepted
% 
%   - value_sizes -
%   * array with the number of elements that each vector of values should
%   have
%
% OUTPUT:
%   - values_array -
%   * Each set of values is returned as an array with a row for each name.
% 
% Original author: Lars D'Hondt
% Original date: 15/April/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

Nfields = length(name_values_cell);
Ncoords = length(valid_names);
Nvalues = length(value_sizes);

% check input size
if mod(Nfields,Nvalues+1)~=0
    error(['Cell array has ' num2str(Nfields) ' inputs, which is not a multiple of ' num2str(Nvalues+1) '.'])
end

% create output
values_array = cell(1,Nvalues);
for k=1:Nvalues
    values_array{k} = nan(value_sizes(k),Ncoords);
end

for ii=1:Nfields/(Nvalues+1)
    names_i = name_values_cell{(Nvalues+1)*(ii-1)+1};
    if ~iscell(names_i)
        names_i = {names_i};
    end
    if strcmp(names_i{1},'all')
        names_i = valid_names;
    end
    for jj=1:length(names_i)
        idx_ij = find(contains(valid_names(:),names_i(jj)));
        if isempty(idx_ij)
            error(['"' names_i{jj} '" is not a valid name.'])
        end
        for k=1:Nvalues
            if isnan(values_array{k}(1,idx_ij))
                if ~isempty(name_values_cell{(Nvalues+1)*(ii-1)+1+k})
                    values_array{k}(:,idx_ij) = vertcat(name_values_cell{(Nvalues+1)*(ii-1)+1+k});
                end
            else
                error(['Multiple values assigned to "' valid_names{idx_ij} '".'])
            end
        end
    end
end

for k=1:Nvalues
    varargout{k} = values_array{k};
end 


