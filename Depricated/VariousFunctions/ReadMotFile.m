function [out] = ReadMotFile(fpath, title)
%SIMM_ReadMotion  Read SIMM motion file.
%
%   Syntax:
%    [header, names, data, fpath] = SIMM_ReadMotion(fpath, title)
%
%   Input:
%    fpath: String containing a path to the file that is to be read. If set
%           to an empty string, the path will be obtained through a file
%           input dialog. Optional, defaults to an empty string.
%    title: Window title for the file input dialog. Only used when fpath is
%           set to an empty string. Optional, defaults to 'Open motion
%           file' if omitted or empty.
%
%   Output:
%    header: String containing the file's header. Will be empty if the user
%            closes the file open dialog box.
%    names:  N-by-1 cell array of strings containing the variable names for
%            which the motion file contains values. Will be empty if the
%            user closes the file open dialog box.
%    data:   M-by-N array containing the trajectory data for all variables
%            in the motion file. Each row corresponds to one time frame;
%            each column corresponds to a cell in names. Will be empty if
%            the user closes the file open dialog box.
%    fpath:  String containing a path to the file that was read. Will be
%            empty if the user closes the file open dialog box.
%
%   Effect: This function will read a SIMM motion file. The header is split
%   off from the data section and returned as an unmodified string; the
%   data section is parsed into a cell array of strings for the variable
%   names and a regular array for the values.
%
%   Dependencies: none
%
%   Known parents: SIMM_ReadMoments.m
%                  OpenSim_ReadSelectSync.m
%                  OpenSim_ReadSOData.m

%Created on 06/05/2008 by Ward Bartels.
%WB, 09/12/2011: Carriage return character is now discarded.
%Stabile, fully functional.


%Set defaults for input variables
if nargin<1, fpath = ''; end
if nargin<2 || isempty(title), title = 'Open motion file'; end

%If no file path was provided, get file from file input dialog
if isempty(fpath)
    [filename, pathname] = uigetfile({'*.mot' 'Motion file (*.mot)'}, title);
    if isequal(filename, 0)
        header = '';
        names = cell(0, 1);
        data = [];
        return
    end
    fpath = fullfile(pathname, filename);
end

%Catch errors while reading string, attempt to close file if one was caught
try
    fid = fopen(fpath, 'rt');
    str = fread(fid, [1 Inf], '*char');
    fclose(fid);
catch err
    if ~isempty(fopen(fid))
        fclose(fid);
    end
    rethrow(err);
end

%Remove comments from string
str = regexprep(str, '/\*.*?\*/', '');

%Split off header
[start, stop] = regexpi(str, '\n\s*endheader\s*\n', 'start', 'end', 'once');
if numel(start)~=1
    error([mfilename ':HeaderEnd'], 'Unique header end not found.');
end
header = str(1:start);
str = str(stop:end);

%Get column names
[names, pos] = textscan(str, '%s', 1, 'whitespace', '\n');
if isempty(names{1})
    error([mfilename ':ColNames'], 'Column name row not found.');
end
names = textscan(names{1}{1}, '%s');
names = names{1};

%Read numerical data
if title == 1
    data = textscan(str(pos:end), '%n');
    data = data{1};
    if mod(numel(data), numel(names)-2)~=0
        error([mfilename ':numdata'], 'incorrect number of data values.');
    end
    data = reshape(data, numel(names)-2, []).';
else
    data = textscan(str(pos:end), '%n');
    data = data{1};
    if mod(numel(data), numel(names))~=0
        error([mfilename ':numdata'], 'incorrect number of data values.');
    end
    data = reshape(data, numel(names), []).';
end

% output
out.header = header;
out.names = names;
out.data = data;
out.fpath = fpath;