function [S,model_info]=preprocessing(S,osim_path)

[]=osim2dll(S,osim_path);

model_info = get_model_info(S,osim_path);

model_info.muscle_info.params = read_and_scale_MTparameters(S,osim_path,model_info);

IOfields = fields(model_info.ExtFunIO);
for i=1:length(IOfields)
    model_info.ExtFunIO.(IOfields{i}) = convert2double(model_info.ExtFunIO.(IOfields{i}));
end
[model_info,S] = GetIndexHelper(S,model_info);

[model_info] = scale_MTparameters(S,model_info);

[model_info] = get_musculoskeletal_geometry_approximation(S,osim_path,model_info);

model_info = update_model_info_muscleProperties(model_info);
