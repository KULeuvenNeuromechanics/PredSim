function [] = write_json(dataStruct, jsonFile)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% tables are unfortunatly not supported. convert these to matrices
% temporarely

% loop over fieldnames to check for tables
fields = fieldnames(dataStruct);
for f = 1:length(fields)
    if istable(dataStruct.(fields{f}))
        % convert table to array
        dataStruct.(fields{f}) = table2array(dataStruct.(fields{f}));
    end
end


% Convert structure to JSON text
jsonText = jsonencode(dataStruct);

% % Optionally pretty-print (for readability)
jsonText = prettyjson(jsonText);  % helper function below

% check if we have to create output folder
[folder, ~, ~] = fileparts(jsonFile);
if ~isfolder(folder)
    mkdir(folder)
end

% Write to file
fid = fopen(jsonFile, 'w');
if fid == -1
    error('Cannot create JSON file.');
end
fwrite(fid, jsonText, 'char');
fclose(fid);


end