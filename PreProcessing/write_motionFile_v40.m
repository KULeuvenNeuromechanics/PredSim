function write_motionFile_v40(q, fname)
% --------------------------------------------------------------------------
% write_motionFile_v40 
%   Write a .mot file compatible with OpenSim 4
% 
% INPUT:
%   q
%   * Structure of data
%  
%   fname
%   * Filename that you want to write, should contain '.mot' (e.g. gait.mot)
% 
% OUTPUT:
%   NA
% 
% Original author: Dhruv Gupta
% Original date: 12/02/21
%
% Last change: Lars D'Hondt
% Change date: 21/02/22
% Changelog:
% If an error occurs while writing, it will now close the file
% --------------------------------------------------------------------------

fid = fopen(fname, 'w');	
if fid == -1								
	error(['unable to open ', fname])		
end

if length(q.labels) ~= size(q.data,2)
	error('Number of labels doesn''t match number of columns')
end

if q.labels{1} ~= 'time'
	error('Expected ''time'' as first column')
end

try
    fprintf(fid, '%s\n', fname);
    fprintf(fid,'version=1\n');
    fprintf(fid, 'nRows=%d\n', size(q.data,1));
    fprintf(fid, 'nColumns=%d\n', size(q.data,2));
    fprintf(fid, 'inDegrees=%s\n', q.inDeg);
    fprintf(fid, 'endheader\n');
    
    for i=1:length(q.labels)
	    fprintf(fid, '%20s\t', q.labels{i});
    end
    fprintf(fid, '\n');
    
    for i=1:size(q.data,1)
	    fprintf(fid, '%20.8f\t', q.data(i,:));
	    fprintf(fid, '\n');
    end
    
catch ME
    fclose(fid);
    error(ME)
end

fclose(fid);
return;
