function [] = AddCasadiPaths()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
name = getenv('COMPUTERNAME');
if strcmp(name,'GBW-D-W2711')
    addpath(genpath('C:\GBW_MyPrograms\casadi-windows-matlabR2016a-v3.5.1'));
elseif strcmp(name,'GBW-L-W0696')
    addpath(genpath('C:\Users\u0088756\Documents\Software\Casadi351'));;
end


end

