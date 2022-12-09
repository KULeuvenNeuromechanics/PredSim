function [ExePath] = InstallOsim2Dll_Exe(ExeDir)
% --------------------------------------------------------------------------
% InstallOsim2Dll_Exe
%   This functions Downloads the compiled opensimAD python utilities to
%   convert .osim model to a .cpp file. The python code can be found in the
%   opensimAD submodule (see opensimAD/osimtocppexe.py) or on github
%   https://github.com/Lars-DHondt-KUL/opensimAD.
%   Previous versions of the compiled utilities can be found on
%   https://www.dropbox.com/sh/05d7tbk9x7dhldq/AAATjYa9QDjdsjAR3QWyBgtta?dl=0 
% 
% INPUT:
%   - ExeDir -
%   * path to folder where you want to install the executable.
%
% OUTPUT:
%   - ExePath -
%   * path to the folder with the executable (This is input setting
%   S.Cpp2Dll.PathCpp2Dll_Exe in run_pred.m)
% 
% Original author: Maarten Afschrift
% Original date: 17/May/2022
%
%   update: added functionality to get new executables if they are outdated
% Last edit by: Lars D'Hondt
% Last edit date: 16/Nov/2022
% --------------------------------------------------------------------------

% link to binaries
Link_ExeZip = 'https://www.dropbox.com/sh/05d7tbk9x7dhldq/AAAUtC7quxwvFtlXFA-R5dqfa/Cpp2Dll_Bin.zip?dl=1'; 
% date of last update of binaries (REQUIRES MANUAL UPDATING)
Date_ExeZip = datetime([2022 12 05],"Format",'uuuu-MM-dd');

% note: I'm aware this is a crude solution. Feel free to suggest a way to
% get matlab to know when a file on dropbox was last updated.

% test if we have to create a new folder
if ~isfolder(ExeDir)
    mkdir(ExeDir)
end

% path to the .zip file
ExePath = fullfile(ExeDir,'Cpp2Dll_Bin');

if isfolder(ExePath) && isfile(fullfile(ExePath,'osimtocppexe.exe')) && isfile(fullfile(ExePath,'genF.exe'))
    dir_osimtocppexe = dir(fullfile(ExePath,'osimtocppexe.exe'));
    is_newer_osim2cpp = datetime(dir_osimtocppexe(1).datenum,'ConvertFrom','datenum') > Date_ExeZip;
    dir_genF = dir(fullfile(ExePath,'genF.exe'));
    is_newer_genF = datetime(dir_genF(1).datenum,'ConvertFrom','datenum') > Date_ExeZip;
    if is_newer_osim2cpp && is_newer_genF
        disp('Program to convert .osim model to .dll file was already downloaded');
        download_bin = 0;
    else
        disp('Updating Osim2Dll program to newer version');
        download_bin = 1;
        delete(fullfile(ExeDir,'Cpp2Dll_Bin.zip'));
        rmdir(fullfile(ExeDir,'Cpp2Dll_Bin'),'s');
    end
else
    disp('Start installation Osim2Dll program');
    download_bin = 1;
end


if download_bin 
    % download the executable from dropbox
    disp('   downloading executable from dropbox');
    websave(fullfile(ExeDir,'Cpp2Dll_Bin.zip'),Link_ExeZip);
    disp('   unzip bin');
    unzip(fullfile(ExeDir,'Cpp2Dll_Bin.zip'),ExeDir);

    % install opensim binaries if needed
    if ~isfolder('C:\OpenSim 4.3\bin') || ~isfile('C:\OpenSim 4.3\bin\OpenSim64.exe')
        disp('   installing opensim binaries in C:\OpenSim 4.3\bin')
        LinkOsimLibraries = 'https://www.dropbox.com/s/x6vvzltb17mtenc/bin.zip?dl=1';        
        websave(fullfile(ExeDir,'OsimBin.zip'),LinkOsimLibraries);
        if ~isfolder('C:\OpenSim 4.3')
            mkdir('C:\OpenSim 4.3');
        end
        unzip(fullfile(ExeDir,'OsimBin.zip'),'C:\OpenSim 4.3');
    end
    disp('   installation Osim2Dll program finished');
end

disp(' ')


end