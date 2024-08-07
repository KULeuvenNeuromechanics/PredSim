clear
close all
clc

%%
casadi_path = casadi.GlobalOptions.getCasadiPath();

computer_options = {'casadi_path',casadi_path};


test1_options = {'max_iter',5, 'compare_result',false};
test2_options = {'compare_printout',false, 'run_as_batch_job',true};


test_models = {'Falisse_et_al_2022', 'DHondt_et_al_2024_4seg'};

%%

for model_name=string(test_models)
    runTestSimulation(model_name, computer_options{:}, test1_options{:});
end
