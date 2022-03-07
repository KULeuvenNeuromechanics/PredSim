function [model_info] = get_musculoskeletal_geometry_approximation(S,osim_path,model_info)

muscle_data = muscle_analysis(S,osim_path,model_info);

[model_info] = PolynomialFit(muscle_data);
