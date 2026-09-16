function velStrings = velocityToString(vel)
% velocityToString Converts velocity values to string format 'XpYms'
%   Input:  vel - numeric scalar or vector (e.g. [0.7 1.1])
%   Output: velStrings - cell array of strings (e.g. {'0p7ms', '1p1ms'})

    % Ensure input is numeric
    if ~isnumeric(vel)
        error('Input must be numeric.');
    end

    % Preallocate cell array for output
    velStrings = cell(size(vel));

    % Loop through each element
    for i = 1:numel(vel)
        % Convert number to string, replace '.' with 'p'
        str = strrep(num2str(vel(i)), '.', 'p');
        % Append 'ms'
        velStrings{i} = [str 'ms'];
    end
end