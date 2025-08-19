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

% loop over inputs
for ii=1:Nfields/(Nvalues+1)
    names_i = name_values_cell{(Nvalues+1)*(ii-1)+1};
    if ~iscell(names_i)
        names_i = {names_i};
    end
    if strcmp(names_i{1},'all')
        names_i = valid_names;
    end
    for jj=1:length(names_i)
        % test that the 1st input is a name
        idx_ij = find(contains(valid_names(:),names_i(jj)));
        if isempty(idx_ij)
            error(['"' names_i{jj} '" is not a valid name.'])
        end
        for i=1:length(idx_ij)
            idx_iji = idx_ij(i);
            % assigns the input values to output arrays for this name
            for k=1:Nvalues
                % test that we have not yet set this array value
                if isFirstEntry(values_array{k}(1,idx_iji))
                    value_k = name_values_cell{(Nvalues+1)*(ii-1)+1+k};
                    if ~isempty(value_k)
                        % if the value is text, convert values_array{k}
                        if isa(value_k,'char')
                            value_k = string(value_k);
                        end
                        if isa(value_k,'string') && ~isa(values_array{k},'string')
                            values_array{k} = string(values_array{k});
                        end
                        values_array{k}(:,idx_iji) = vertcat(value_k(:));
                    end
                else
                    error(['Multiple values assigned to "' valid_names{idx_iji} '".'])
                end
            end
        end
    end
end

for k=1:Nvalues
    varargout{k} = values_array{k};
end 


end
% helper function
function [has_no_entry] = isFirstEntry(test_field)
    has_no_entry = 0;
    if isa(test_field,'double') && isnan(test_field)
        has_no_entry = 1;
    elseif isa(test_field,'string') && ismissing(test_field)
        has_no_entry = 1;
    end
end