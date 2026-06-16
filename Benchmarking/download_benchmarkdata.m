function [outfolder] = download_benchmarkdata(varargin)
%download_predsimdata Downloads data from https://github.com/MaartenAfschrift/predsim_benchmark_data
% and saves it in PredSim/Benchmarking/Data

if ~isempty(varargin)
    bool_overwrite = varargin{1};
else
    bool_overwrite = true;
end

% get current path
fullPath = mfilename('fullpath');
folderPath = fileparts(fullPath);
outfolder = fullfile(folderPath,'data');
if ~isfolder(outfolder)
    mkdir(outfolder);
end

if bool_overwrite || ~exist(fullfile(outfolder,'Readme.md'),'file')

    % download zip file with data
    websave(fullfile(outfolder,'temp.zip'),...
        ['https://github.com/MaartenAfschrift/predsim_benchmark_data/', ....
        'archive/refs/heads/master.zip']);

    % unpack zip file
    unzip(fullfile(outfolder,'temp.zip'),outfolder);
    delete(fullfile(outfolder,'temp.zip'));

    % copy past some things
    src = fullfile(outfolder, 'predsim_benchmark_data-master', '*');  % note the *
    copyfile(src, outfolder);
    rmdir(fullfile(outfolder,'predsim_benchmark_data-master'), 's');
end








end