function [S] = adaptExport3DSegmentOrigins(S,osim_path)
% --------------------------------------------------------------------------
% adaptExport3DSegmentOrigins
%   Remove the entries of S.Cpp2Dll.export3DSegmentOrigins that do not
%   exist in the osim model.
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%
% OUTPUT:
%   - S -
%   * setting structure S
%
% 
% Original author: Lars D'Hondt
% Original date: 6/Dec/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

segments = S.Cpp2Dll.export3DSegmentOrigins;

if isempty(segments)
    return
end

import org.opensim.modeling.*;
model = Model(osim_path);
bodyset = model.getBodySet();

for i=1:length(segments)
    try
        body = bodyset.get(segments{i});
    catch
        disp(['   Cannot find body names "' segments{i} '" in osim model. ',...
            'Removing name from S.Cpp2Dll.export3DSegmentOrigins'])
        segments{i} = 'invalid';
    end
end

segments = segments(~strcmp(segments,'invalid'));

S.Cpp2Dll.export3DSegmentOrigins = segments;

end