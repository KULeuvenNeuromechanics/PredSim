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