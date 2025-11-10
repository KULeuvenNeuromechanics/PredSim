%% Test benchmark plot function

% path with benchmark results
resfolder = ['C:\Users\mat950\OneDrive - Vrije Universiteit Amsterdam\' ...
    'Onderzoek\SimResults\Benchmark_Falisse2022_v2'];

% I have to adapt a file slightly due to a previous error
benchmark_filepath = fullfile(resfolder,"benchmark_settings.mat");
load(benchmark_filepath,'S','osim_path','S_benchmark');

S_benchmark.koelewijn.names = {'Falisse_et_al_2022_down_8_speed_80','Falisse_et_al_2022_speed_130',...
    'Falisse_et_al_2022_up_8_speed_80','Falisse_et_al_2022_up_8_speed_130',...
    'Falisse_et_al_2022_speed_80','Falisse_et_al_2022_speed_130'};
save(benchmark_filepath,'S','osim_path','S_benchmark');

% default function to evaluate benchmarking results
benchmark_results(resfolder,'BoolPlot',true);


