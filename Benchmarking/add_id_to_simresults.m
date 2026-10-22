function [] = add_id_to_simresults(datapath, id)
%add_id_to_simresults Simple function to add identifier of a simulation to 
% the simulations results. This was mainly used to update old simulations
% so that we don't have to run the whole thing again to simply add an
% identifier to the simulation resultsfile


% get all matfiles
matfiles = dir(fullfile(datapath,'*.mat'));
for ifile = 1:length(matfiles)
    matfile_sel = fullfile(matfiles(ifile).folder, matfiles(ifile).name);
    if ~ismember('R',who('-file',matfile_sel))
        continue
    end
    data = load(matfile_sel,'R');
    R = data.R;
    if isfield(R.S.misc,'benchmark_id') && ~isempty(R.S.misc.benchmark_id)
        if ~strcmp(R.S.misc.benchmark_id,id)
            warning('PredSim:BenchmarkIdConflict',...
                'Keeping existing ID %s in %s; requested ID is %s. Verify the original simulation condition before migrating it.',...
                R.S.misc.benchmark_id,matfile_sel,id);
        end
        continue
    end
    % add id to settings
    R.S.misc.benchmark_id = id;
    save(matfile_sel,'R','-append');
end



end
