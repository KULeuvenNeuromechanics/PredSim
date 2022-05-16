function [S] = osim2dll(S,osim_path)
% --------------------------------------------------------------------------
% osim2dll
%   Matlab shell for the python workflow
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
%
% OUTPUT:
%   - S -
%   * setting structure S
% 
% Original author: Lars D'Hondt
% Original date: 07/March/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


% extract name of opensim model file
[~,osim_file_name,~] = fileparts(osim_path);

% external function name
S.misc.external_function = ['F_' osim_file_name '.dll'];

if ~isfile(fullfile(S.misc.subject_path,S.misc.external_function))
    pathOpenSimModel = osim_path;
    outputFilename = ['F_' osim_file_name];
    compiler = "Visual Studio 15 2017 Win64";
    export3DSegmentOrigins = {'calcn_r', 'calcn_l', 'femur_r', 'femur_l', 'hand_r',...
                          'hand_l', 'tibia_r', 'tibia_l', 'toes_r', 'toes_l'};



end


