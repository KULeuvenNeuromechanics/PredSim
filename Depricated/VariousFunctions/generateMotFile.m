% Generate a motion *.mot file readable by OpenSim
% 
% Tim Dorn
% November 2008
% 
% --------------------------------------------------------------------
% Usage: generateMotFile(dataMatrix, colnames, filename)
% --------------------------------------------------------------------
% 
% Inputs:   dataMatrix = data matrix to write to file
%                       (first column should be time)
%           colnames = cell array of column name strings
%           filename = string containing the output filename (must include extension)
% 
% Outputs:  output motion file
% 
% 
% Notes:    Number of data columns must match the number of column names or
%           an exception will be thrown.
% 
% --------------------------------------------------------------------
% 
% Copyright (c)  2008 Tim Dorn
% Use of the GaitExtract Toolbox is permitted provided that the following
% conditions are met:
% 	1. The software is not distributed or redistributed.  Software distribution is allowed 
%     only through https://simtk.org/home/c3dtoolbox.
% 	2. Use of the GaitExtract Toolbox software must be acknowledged in all publications,
%      presentations, or documents describing work in which the GaitExtract Toolbox was used.
% 	3. Credits to developers may not be removed from source files
% 	4. Modifications of source code must retain the above copyright notice, this list of
%     conditions and the following disclaimer. 
% 
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
%  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
%  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
%  SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
%  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
%  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
%  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%  OR BUSINESS INTERRUPTION) OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
%  WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% --------------------------------------------------------------------

function generateMotFile(dataMatrix, colnames, filename)

[datarows, datacols] = size(dataMatrix);
time = dataMatrix(:,1);
range = [time(1), time(end)];

if length(colnames) ~= datacols
    error('Number of column names (%d) do not match the number of columns in the data (%d)\n', ...
        length(colnames), datacols);
end


% MOT File Header
% ---------------

fid = fopen(filename, 'w');
if fid < 0
    fprintf('\nERROR: %s could not be opened for writing...\n\n', filename);
    return
end

fprintf(fid, '%s\nnRows=%d\nnColumns=%d\n\n', filename, datarows, datacols);
fprintf(fid, 'name %s\ndatacolumns %d\ndatarows %d\nrange %f %f\nendheader\n', ...
    filename, datacols, datarows, range(1), range(2));


% MOT File Body
% -------------
cols = [];
for i = 1:datacols,
    if i == 1
        cols = [cols, colnames{i}];
    else
        cols = [cols, sprintf('\t%s', colnames{i})]; 
    end
end
cols = [cols, '\n'];
fprintf(fid, cols);

for i = 1:datarows,
    fprintf(fid, '%20.10f\t', dataMatrix(i,:));
    fprintf(fid, '\n');
end

fclose(fid);
%fprintf('Saved motion file: %s\n', filename);

