function [St] = convert2double(St)
% --------------------------------------------------------------------------
% convert2double
%   This function converts data in each field in a struct to type double.
%   
% INPUT:
%   - St -
%   * struct with fields that contain int32 values
%
% OUTPUT:
%   - St -
%   * struct with fields that contain double values
% 
% Original author: Dhruv Gupta
% Original date: 4/Feb/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

St_fields = fields(St);
for i=1:length(St_fields)
    St.(St_fields{i}) = double(St.(St_fields{i}));
end
    
