function exists = benchmark_result_exists(folder)
% A settings file or failed/incomplete solve must not suppress a retry.
exists = false;
files = dir(fullfile(folder,'*.mat'));
for k = 1:numel(files)
    R = load_sim_file(fullfile(files(k).folder,files(k).name));
    if ~isempty(R) && all(isfield(R,{'kinematics','spatiotemp','metabolics'}))
        exists = true;
        return
    end
end
end
