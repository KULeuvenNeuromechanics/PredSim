function [data] = read_json(jsonFile)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% Read the JSON file as text
fid = fopen(jsonFile, 'r');
if fid == -1
    error('Cannot open JSON file.');
end
rawText = fread(fid, '*char')';
fclose(fid);

% Decode JSON text into a MATLAB structure
data = jsondecode(rawText);

% convert ik and id to a table if possible
% loop over fieldnames to check for tables
fields = fieldnames(data);
for f = 1:length(fields)
    if strcmp(fields{f},'ik') || strcmp(fields{f},'ik_std')
        if ~isempty(data.(fields{f}))
            % convert array to table
            data.(fields{f}) = array2table(data.(fields{f}),...
                'VariableNames',data.ik_header);
        end
    elseif strcmp(fields{f},'id') || strcmp(fields{f},'id_std')
        if ~isempty(data.(fields{f}))
            % convert array to table
            data.(fields{f}) = array2table(data.(fields{f}),...
                'VariableNames',data.ik_header);
        end
    end
end





end