function [] = run_pred_sim(S,osim_path)

S = getDefaultSettings(S);
S = updateS_passiveJointTorqueProperties(S);

[S,model_info] = preprocessing(S,osim_path);

[f_casadi] = createCasadiFunctions(model_info);

OCP_formulation(S,model_info,f_casadi);

post_processing(S,model_info,f_casadi);
