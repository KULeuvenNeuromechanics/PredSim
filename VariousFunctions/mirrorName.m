function [varargout] = mirrorName(original_name)
% --------------------------------------------------------------------------
% mirrorName
%   This functions creates the name of the mirrored element (joint, muscle,
%   etc.) based on the input. 
%   Examples:
%       ans1 = mirrorName('Soleus_r') 
%           ans1 = 'Soleus_l'
%
%       [ans1, ans2] = mirrorName('R_knee')
%           ans1 = 'L_knee'
%           ans2 = 'R_knee' 
%
% INPUT:
%   - original_name -
%   * name of the left or right side. Side should be denoted by 'l_*' or
%   '*_R'. Not case sensitive. (char or string)
%
% OUTPUT:
%   - mirrored_name (if 1 output argument) -
%   * name of the right or left side, based on given input (char)
%
%   - [left_name, right_name] (if 2 output arguments) -
%   * name of the left and right side, based on given input (char)
%
% 
% Original author: Lars D'Hondt
% Original date: 17/April/2023
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

mirrored_name = [];
left_name = [];
right_name = [];

original_name = char(original_name);
suffix = original_name(end-1:end);

switch suffix
    case '_r'
        right_name = original_name;
        left_name = [original_name(1:end-2) '_l'];
        mirrored_name = left_name;

    case '_l'
        left_name = original_name;
        right_name = [original_name(1:end-2) '_r'];
        mirrored_name = right_name;

    case '_R'
        right_name = original_name;
        left_name = [original_name(1:end-2) '_L'];
        mirrored_name = left_name;

    case '_L'
        left_name = original_name;
        right_name = [original_name(1:end-2) '_R'];
        mirrored_name = right_name;

    otherwise
        suffix = -1;
end

if suffix == -1

    prefix = original_name(1:2);

    switch prefix
        case 'r_'
            right_name = original_name;
            left_name = ['l_' original_name(3:end)];
            mirrored_name = left_name;

        case 'l_'
            left_name = original_name;
            right_name = ['r_' original_name(3:end)];
            mirrored_name = right_name;

        case 'R_'
            right_name = original_name;
            left_name = ['L_' original_name(3:end)];
            mirrored_name = left_name;

        case 'L_'
            left_name = original_name;
            right_name = ['R_' original_name(3:end)];
            mirrored_name = right_name;

        otherwise
            % not detected as left or right side
            prefix = -1;
            left_name = original_name;
            right_name = original_name;
            mirrored_name = original_name;
    end
end

% return output
if nargout == 1
    varargout{1} = mirrored_name;

elseif nargout == 2
    varargout{1} = left_name;
    varargout{2} = right_name;
end


end % end of function