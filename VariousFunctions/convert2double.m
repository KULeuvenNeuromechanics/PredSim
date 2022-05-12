function [St1] = convert2double(St)
St_fields = fields(St);
for i=1:length(St_fields)
    St1.(St_fields{i}) = double(St.(St_fields{i}));
end
    
