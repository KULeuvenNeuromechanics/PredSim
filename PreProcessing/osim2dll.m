function [S] = osim2dll(S,osim_path)
%
% Matlab shell for the python workflow
%
% 
% Author: Lars D'Hondt
%
% Date: 07/March/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract name of opensim model file
[~,osim_file_name,~] = fileparts(osim_path);

% external function name
S.misc.external_function = ['F_' osim_file_name '.dll'];




