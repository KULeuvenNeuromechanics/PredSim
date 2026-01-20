function [var_is_env] = isenv(varname)
% isenv
%   Determine if environment variable exists.
%   To be used with version before MATLAB R2022b
%   https://mathworks.com/help/matlab/ref/isenv.html
varname = string(varname);
var_is_env = false(size(varname));
for i=1:size(varname,1)
    for j=1:size(varname,2)
        var_is_env(i,j) = ~isempty(getenv(varname(i,j)));
    end
end

end