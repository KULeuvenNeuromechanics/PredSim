function specific_tension = getSpecificTensions(S, muscleNames)
% --------------------------------------------------------------------------
% getSpecificTensions
%   Returns the specific tension value for each muscle.
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - muscleNames -
%   * Cell array of muscle names
%
% OUTPUT:
%   - specific_tension -
%   * Array with specific tension of each muscle. Default is 0.70
% 
% Original author: Lars D'Hondt
% Original date: 16 September 2025
% --------------------------------------------------------------------------

specific_tension = 0.70*ones(length(muscleNames),1);

if exist(S.subject.default_specific_tension,'file')
    default_sigma = readtable(S.subject.default_specific_tension);

    for i=1:length(muscleNames)
        default_sigma_i = default_sigma(strcmp(default_sigma.name, muscleNames{i}));
        if ~isempty(default_sigma_i)
            specific_tension(i) = default_sigma_i;
        end
    end

elseif isnumeric(S.subject.default_specific_tension) && S.subject.default_specific_tension > 0
    specific_tension(:) = S.subject.default_specific_tension;

end

specific_tension = unpack_name_value_combinations(...
    S.subject.set_specific_tension_selected_muscles, muscleNames, 1, specific_tension);


end    
