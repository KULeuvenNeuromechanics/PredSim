function [ExePath] = InstallOsim2Dll_Exe(ExeDir)
% --------------------------------------------------------------------------
% InstallOsim2Dll_Exe
%   This functions Downloads the compiled opensimAD python utilities to
%   convert .osim model to a .cpp file. The python code can be found in the
%   opensimAD submodule (see opensimAD/osimtocppexe.py) or on github
%   https://github.com/Lars-DHondt-KUL/opensimAD
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
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


% test if we have to create a new folder
if ~isfolder(ExeDir)
    mkdir(ExeDir)
end

% path to the .zip file
ExePath = fullfile(ExeDir,'Cpp2Dll_Bin');

if isfolder(ExePath) || isfile(fullfile(ExePath,'osimtocppexe.exe'))
    disp('Program to convert .osim model to .dll file was already downloaded');
else
    % download the executable from google drive
    disp('Start installation osim2Dll program');
    disp('   downloading executable from dropbox');
    Link_ExeZip = 'https://www.dropbox.com/s/43iyndhj2vxr20v/Cpp2Dll_Bin.zip?dl=1';
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
    disp('   installation osim2Dll program finished');
end




end