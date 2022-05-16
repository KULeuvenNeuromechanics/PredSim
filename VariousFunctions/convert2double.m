function [St] = convert2double(St)
% --------------------------------------------------------------------------
% convert2double
%   (Explanation)
%   
% INPUT:
%   - St -
%   * 
%
% OUTPUT:
%   - St -
%   * 
% 
% Original author: 
% Original date: 
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

St_fields = fields(St);
for i=1:length(St_fields)
    St.(St_fields{i}) = double(St.(St_fields{i}));
end
    
