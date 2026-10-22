function test_benchmark_plot(resfolder)
% Plot a benchmark run; default to the accompanying Falisse example output.
if nargin == 0
    repo = fileparts(fileparts(mfilename('fullpath')));
    resfolder = fullfile(repo,'Results','Benchmark_Falisse2022');
end
add_benchmarkdata_to_simresults(resfolder,'BoolPlot',true);
end
