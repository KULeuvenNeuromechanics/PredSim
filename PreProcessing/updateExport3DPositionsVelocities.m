function [S] = updateExport3DPositionsVelocities(S,osim_path)
% --------------------------------------------------------------------------
% updateExport3DPositionsVelocities
%   Updates the entries of S.OpenSimADOptions.export3DPositions and 
%   S.OpenSimADOptions.export3DVelocities. 
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
% Original date: 9/May/2023s
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import org.opensim.modeling.*;

%% Fill in blank inputs in points definitions

for i=1:length(S.bounds.points)
    % body needs to be specified
    if ~isfield(S.bounds.points(i),'body') || isempty(S.bounds.points(i).body)
        errmsg = ['Please enter the name of a body in your osim model in "S.bounds.points(',...
            num2str(i) ').body".'];
        error(errmsg);
    end

    % default position is origin
    if ~isfield(S.bounds.points(i),'point_in_body') || isempty(S.bounds.points(i).point_in_body)
        S.bounds.points(i).point_in_body = [0,0,0];
    end

    % default name is name of body
    if ~isfield(S.bounds.points(i),'name') || isempty(S.bounds.points(i).name)
        S.bounds.points(i).name = S.bounds.points(i).body;
    end

    pointNames{i} = S.bounds.points(i).name;
end

%% Fill in blank inputs in distanceConstraints

for i=1:length(S.bounds.distanceConstraints)
    % points need to be specified
    if ~isfield(S.bounds.distanceConstraints(i),'point1')...
            || isempty(S.bounds.distanceConstraints(i).point1)
        error(['Please enter the name of a point defined in "S.bounds.points", ',...
            'in "S.bounds.distanceConstraints(', num2str(i) ').point1".'])
    end
    if ~isfield(S.bounds.distanceConstraints(i),'point2')...
            || isempty(S.bounds.distanceConstraints(i).point2)
        error(['Please enter the name of a point defined in "S.bounds.points", ',...
            'in "S.bounds.distanceConstraints(', num2str(i) ').point2".'])
    end

    % if point was not defined, assume the name refers to the body and the
    % point lies in the body origin
    if ~any(strcmp(pointNames,S.bounds.distanceConstraints(i).point1))
        S.bounds.points(end+1).body = S.bounds.distanceConstraints(i).point1;
        S.bounds.points(end).point_in_body = [0,0,0];
        S.bounds.points(end).name = S.bounds.distanceConstraints(i).point1;
    end
    if ~any(strcmp(pointNames,S.bounds.distanceConstraints(i).point2))
        S.bounds.points(end+1).body = S.bounds.distanceConstraints(i).point2;
        S.bounds.points(end).point_in_body = [0,0,0];
        S.bounds.points(end+1).name = S.bounds.distanceConstraints(i).point2;
    end

    % by default, take 3D distance
    if ~isfield(S.bounds.distanceConstraints(i),'direction')...
            || isempty(S.bounds.distanceConstraints(i).direction)
        S.bounds.distanceConstraints(i).direction = 'xyz';
    end

    % constraint needs at least 1 bound
    hasLowerBound = isfield(S.bounds.distanceConstraints(i),'lower_bound')...
            && ~isempty(S.bounds.distanceConstraints(i).lower_bound);
    hasUpperBound = isfield(S.bounds.distanceConstraints(i),'upper_bound')...
            && ~isempty(S.bounds.distanceConstraints(i).upper_bound);
    if ~hasLowerBound && ~hasUpperBound
        error(['The distance constraint between "' S.bounds.distanceConstraints(i).point1, ...
            '" and "' S.bounds.distanceConstraints(i).point2 '" ("S.bounds.distanceConstraints(',...
            num2str(i) ')") has no value for upper or lower bound.'])
    end


end


%% Add points for constraints to export3DPositions

S.OpenSimADOptions.export3DPositions = [S.OpenSimADOptions.export3DPositions(:), S.bounds.points(:)];


%% Test that the bodies exist in the osim model

fields = ["export3DPositions","export3DVelocities"];

model = Model(osim_path);
bodyset = model.getBodySet();

for j=fields
    segments = S.OpenSimADOptions.(j);
    
    if isempty(segments)
        continue
    end

    validNames = zeros(length(segments,1));
    
    for i=1:length(segments)
        try
            body = bodyset.get(segments(i).body);
            validNames(i) = 1;
        catch
            warning(['   Cannot find body name "' segments(i).body '" in osim model. ',...
                'Removing from S.OpenSimADOptions.' char(j) '.'])
        end
    end
    
    segments = segments(validNames);
    
    S.OpenSimADOptions.(j) = segments;
end

%% Remove constraints that rely on undefined points

validConstraints = zeros(length(S.bounds.distanceConstraints),1);
for i=1:length(S.OpenSimADOptions.export3DPositions)
    pointNames(i) = S.OpenSimADOptions.export3DPositions(i).name;
end
for i=1:length(S.bounds.distanceConstraints)

    if ~any(strcmp(pointNames,S.bounds.distanceConstraints(i).point1))
        validConstraints(i) = 0;
        warning(['   Cannot find valid point "' S.bounds.distanceConstraints(i).point1,...
            '". Removing entry from S.bounds.distanceConstraints.'])

    elseif ~any(strcmp(pointNames,S.bounds.distanceConstraints(i).point2))
        validConstraints(i) = 0;
        warning(['   Cannot find valid point "' S.bounds.distanceConstraints(i).point2,...
            '". Removing entry from S.bounds.distanceConstraints.'])
    end


end


%% Add field with indices of direction of distance
for i=1:length(S.bounds.distanceConstraints)

    dir_i = char(S.bounds.distanceConstraints(i).direction);

    if length(dir_i)<=3
        idx_dir = [0,0,0];
        if contains(dir_i,'x')
            idx_dir(1) = 1;
        end
        if contains(dir_i,'y')
            idx_dir(2) = 1;
        end
        if contains(dir_i,'z')
            idx_dir(3) = 1;
        end
    else

        switch lower(dir_i)
            case 'sagittal'
                idx_dir = [1,1,0];
            case {'coronal','frontal'}
                idx_dir = [0,1,1];
            case {'transverse'}
                idx_dir = [1,0,1];

            otherwise
                error(['Direction of distance constraint is invalid. Please set',...
                    ' "S.bounds.distanceConstraints(' num2str(i) ').direction" to any',...
                    ' combination of x, y and z, or an anatomical plane.'])
        end
    end

    S.bounds.distanceConstraints(i).directionVectorIdx = idx_dir;

end













end