function [] = DispHeader(header)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(header)
    disp([num2str(i) ' - ' header{i}]);
end

