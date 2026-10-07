function [R, benchmark, model_info] = load_sim_file(sim_res_folder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Accept an exact file when callers have already enumerated results.
if isfile(sim_res_folder)
    mat_files = dir(sim_res_folder);
else
    mat_files = dir(fullfile(sim_res_folder,'*.mat'));
end
R = [];
benchmark = [];
model_info = [];
% Ignore settings and intermediate MAT files when selecting a result.
is_result = false(size(mat_files));
for k = 1:numel(mat_files)
    is_result(k) = ismember('R',who('-file',fullfile(mat_files(k).folder,mat_files(k).name)));
end
mat_files = mat_files(is_result);
if ~isempty(mat_files)
    if length(mat_files) > 1
        disp('warning mutiple mat files in folder')
        disp(sim_res_folder);
        disp(['assumes that file ' mat_files(1).name, ...
            'contains the simulation results'])
    end
    sim_res_file = fullfile(mat_files(1).folder, mat_files(1).name);
    vars = who('-file',sim_res_file);
    vars = intersect(vars,{'R','benchmark','model_info','stats'});
    data = load(sim_res_file,vars{:});
    if isfield(data,'stats') && isfield(data.stats,'success') && ~data.stats.success
        warning('PredSim:UnconvergedBenchmarkResult','Skipping unconverged result %s.',sim_res_file);
        return
    end
    R = data.R;
    if isfield(data,'benchmark')
        benchmark = data.benchmark;
    end
    if isfield(data,'model_info')
        model_info = data.model_info;
    end
else
    R = [];
    benchmark = [];
    model_info = [];
end

end
