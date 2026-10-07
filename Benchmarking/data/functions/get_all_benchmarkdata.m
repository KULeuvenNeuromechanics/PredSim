function [data,studyList, identifierList, slopeList, speedList] = get_all_benchmarkdata(varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(varargin)
    datafolder = varargin{1};
else
    % Get the full path of the currently running function
    fullPath = mfilename('fullpath');
    % Extract the folder path
    folderPath = fileparts(fullPath);
    datafolder = fileparts(folderPath);
end


% find all json files
files = dir(fullfile(datafolder, '**', '*.json'));

% read all data
data = cell(0);
for i =1:length(files)
    data{i} = read_json(fullfile(files(i).folder, files(i).name));
end

% get some lists that make it easier for indexing
studyList = cellfun(@(x) x.study, data, 'UniformOutput', false);
identifierList =  cellfun(@(x) x.identifier, data, 'UniformOutput', false);
slopeList = cellfun(@(x) x.slope, data, 'UniformOutput', false);
speedList = cellfun(@(x) x.speed, data, 'UniformOutput', false);

end