function write_motionFile_v40(q, fname)
% --------------------------------------------------------------------------
% write_motionFile_v40
%   Function to save data to a motion (.mot) file.
%   
% INPUT:
%   - q -
%   * q.labels: cell array with label for each data column
%   * q.data: array with data where each row is a timeframe
%   * q.inDeg: are the rotations in degrees? 'yes' or 'no'
%
%   - fname -
%   * full path and name of the motion file to save
% 
% Original author: Dhruv Gupta
% Original date: 12/02/21
%
% Last edit by: 
% Last edit date: 
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
    error(ME.message)
end

fclose(fid);
return;