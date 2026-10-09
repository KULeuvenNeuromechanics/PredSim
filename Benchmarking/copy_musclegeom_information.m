function [] = copy_musclegeom_information(or_model_path,output_model_path,S)
%copy_musclegeom_information copies muscle geometry information file from
%or_model_path to output_model_path
%   Detailed explanation goes here

% Match getDefaultSettings rather than assuming the default polynomial order.
lower_order = 3;
upper_order = 9;
if nargin >= 3 && isfield(S,'misc') && isfield(S.misc,'poly_order')
    if isfield(S.misc.poly_order,'lower')
        lower_order = S.misc.poly_order.lower;
    end
    if isfield(S.misc.poly_order,'upper')
        upper_order = S.misc.poly_order.upper;
    end
end
suffix = sprintf('_f_lMT_vMT_dM_poly_%g_%g.casadi',lower_order,upper_order);

% folder original model
[folder_or,modelname,~] = fileparts(or_model_path);
% Preprocessing stores caches by subject, even for an external input model.
if nargin >= 3 && isfield(S,'misc') && isfield(S.misc,'main_path') && ...
        isfield(S,'subject') && isfield(S.subject,'name')
    folder_or = fullfile(S.misc.main_path,'Subjects',S.subject.name);
end
geom_file = fullfile(folder_or,[strrep(modelname,' ','_') suffix]);

% folder output
[folder_output,modelname,~] = fileparts(output_model_path);
geom_file_out = fullfile(folder_output,[strrep(modelname,' ','_') suffix]);


% copy file
copyfile(geom_file,geom_file_out);



end
