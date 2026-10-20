function q = read_motionFile_v40(fname)
% --------------------------------------------------------------------------
% read_motionFile_v40
%   This function reads a file in the format of a SIMM motion file
%   and returns a data structure
%   
% INPUT:
%   - fname -
%   * fname is the name of the ascii datafile to be read 
%           ('character array')
%           inDeg: are the values in degrees
%
% OUTPUT:
%   - q -
%   * q returns a structure with the following format:
%				q.labels 	= array of column labels
%				q.data 		= matrix of data
%				q.inDeg 	= values in degrees or radians
%				q.nr 		= number of matrix rows
%				q.nc 		= number of matrix columns
% 
% Original author: ASA 
% Original date: 12/03
%
% Edit by: Eran Guendelman
% Edit date: 09/06
%
% Last edit by: Dhruv Gupta
% Last edit date: 12/02/21
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

fid = fopen(fname, 'r');	
if fid == -1								
	error(['unable to open ', fname])		
end
% Process the file header;
% store # data rows, # data columns.
q.nr = 0; % Added to ensure that the q structures from reading a motion file
q.nc = 0; % are always the same, even if nr and nc are different orders in file.
nextline = fgetl(fid);	
while ~strncmpi(nextline, 'endheader', length('endheader'))
	if strncmpi(nextline, 'nRows', length('nRows'))
		q.nr = str2num(nextline(findstr(nextline, '=')+1 : length(nextline)));
	elseif strncmpi(nextline, 'nColumns', length('nColumns'))
		q.nc = str2num(nextline(findstr(nextline, '=')+1 : length(nextline)));
	elseif strncmpi(nextline, 'inDegrees', length('inDegrees'))
		q.inDeg = nextline(findstr(nextline, '=')+1 : length(nextline));
	end
	nextline = fgetl(fid);
end
% Process the column labels.
nextline = fgetl(fid);
if (all(isspace(nextline))) % Blank line, so the next one must be the one containing the column labels
	nextline = fgetl(fid);
end
q.labels = cell(1, q.nc);
for j = 1:q.nc
	[q.labels{j}, nextline] = strtok(nextline);
end
% Process the data.
% Note:  transpose is needed since fscanf fills columns before rows.
q.data = fscanf(fid, '%f', [q.nc, q.nr])';
fclose(fid);
return;