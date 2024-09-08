% --------------------------------------------------------------------------
% test_PredSim
%   Run tests to ensure PredSim works correctly.
%
%   See also runTestSimulation
%
% Original author: Lars D'Hondt
% Original date: 30 August 2024
% --------------------------------------------------------------------------

clear
close all
clc

%%

% test 1: run simulation for 5 iterations, and compare solver printout
test1_options = {'max_iter',5, 'compare_result',false};

% test 2: run simulation untill convergence, and compare results
test2_options = {'compare_printout',false, 'run_as_batch_job',true};

% models used for reference simulations
test_models = {'Falisse_et_al_2022', 'DHondt_et_al_2024_4seg'};

%%

for model_name=string(test_models)
    runTestSimulation(model_name, test1_options{:});
end
