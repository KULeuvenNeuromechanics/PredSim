function [varargout] = runPredSim(S, osim_path)
% --------------------------------------------------------------------------
% runPredSim
%   Function to run PredSim.
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%
% OUTPUT:
%   - savename (optional) -
%   * results of the simulation will be saved in a file with this name
% 
% Original author: Lars D'Hondt
% Original date: 14/August/2024
% --------------------------------------------------------------------------


if S.solver.run_as_batch_job
    % add to batch
    add_pred_sim_to_batch(S,osim_path)
else
    % run
    [savename] = run_pred_sim(S,osim_path);
end

if nargout == 1
    varargout{1} = savename;
end

end
