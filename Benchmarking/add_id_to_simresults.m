function [] = add_id_to_simresults(datapath, id)
%add_id_to_simresults Simple function to add identifier of a simulation to 
% the simulations results. This was mainly used to update old simulations
% so that we don't have to run the whole thing again to simply add an
% identifier to the simulation resultsfile


% get all matfiles
matfiles = dir(fullfile(datapath,'*.mat'));
for ifile = 1:length(matfiles)
    matfile_sel = fullfile(matfiles(ifile).folder, matfiles(ifile).name);
    load(matfile_sel,'R');
    % add id to settings
    R.S.misc.benchmark_id = id;
    save(matfile_sel,'R','-append');
end



end