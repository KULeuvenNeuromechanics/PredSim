function [MSVS_compiler] = findVisualStudioInstallation()
% --------------------------------------------------------------------------
% findVisualStudioInstallation
%   This functions looks for a Visual Studio installation to fill in
%   S.Cpp2Dll.compiler. If nothing is found, thrown an error to prompt the
%   user to specify this setting manually.
% 
% INPUT:
%   - (none) -
%
% OUTPUT:
%   - MSVS_compiler -
%   * setting structure S
% 
% Original author: Lars D'Hondt
% Original date: 28/November/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

MSVS_years = {'2026','2022','2019','2017','2015'};
MSVS_versions = {'18','17','16','15','14'};
MSVS_bit = {'','','',' Win64',' Win64'}; % 2019 and newer versions assume 64 bit
MSVS_compiler = [];

expected_path = 'c:/Program Files*/Microsoft Visual Studio/';

if ~isempty(dir(expected_path))
    for i=1:length(MSVS_years)
        if ~isempty(dir(fullfile(expected_path,MSVS_years{i})))
            MSVS_compiler = ['Visual Studio ',MSVS_versions{i},' ',MSVS_years{i}, MSVS_bit{i}];
            break
        end
        if ~isempty(dir(fullfile(expected_path,MSVS_versions{i})))
            MSVS_compiler = ['Visual Studio ',MSVS_versions{i},' ',MSVS_years{i}, MSVS_bit{i}];
            break
        end
    end
end

if isempty(MSVS_compiler)
    error(['Could not detect Visual Studio in c:/Program Files or c:/Program Files (x86). ',...
        'Please set S.OpenSimADOptions.compiler to Visual Studio 14 2015 Win64, Visual Studio 15 2017 Win64,',...
        'Visual Studio 16 2019, Visual Studio 17 2022, or Visual Studio 18 2026 based on your installed version']);
end
