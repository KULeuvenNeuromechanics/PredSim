clear
clc
[pathTests,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathTests);
addpath([pathRepo '\PreProcessing'])
addpath([pathRepo '\CasadiFunctions'])


muscle_data = muscleAnalysisAPI(S,osim_path,model_info,50);

[model_info] = PolynomialFit(S,muscle_data,model_info);

f_lMT_vMT_dM_1 = createCasadi_MSKGeometry(S,model_info);


