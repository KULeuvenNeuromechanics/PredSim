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
%
% --------------------------------------------------------------------------
% This file is part of PredSim.
% 
% PredSim: A Framework for Rapid Predictive Simulations of Locomotion
% Copyright (c) 2026 KU Leuven
% 
% PredSim is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Affero General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% 
% PredSim is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public 
% License for more details.
% 
% You should have received a copy of the GNU Affero General Public License 
% along with PredSim. If not, see <https://www.gnu.org/licenses/>.
% --------------------------------------------------------------------------


% Make sure casadi path is set up correctly. This needs to happen before
% adding a simulation to the batch
if ~isfield(S.solver,'CasADi_path')
    try
        S.solver.CasADi_path = casadi.GlobalOptions.getCasadiPath();
    catch
        error("Please add CasADi to the matlab search path, or pass the path " + ...
            "to your CasADi installation (top folder) to S.solver.CasADi_path.")
    end
elseif ~isempty(S.solver.CasADi_path) && ~isfolder(S.solver.CasADi_path)
    error("Unable to find the path assigned to S.solver.CasADi_path:" + ...
        " \n\t%s",S.solver.CasADi_path)
end

% Disable OpenSim log files
org.opensim.modeling.Logger.setLevelString("off");

if S.solver.run_as_batch_job
    % add to batch
    [savename] = add_pred_sim_to_batch(S,osim_path);
else
    % run
    [savename] = run_pred_sim(S,osim_path);
end

if nargout == 1
    varargout{1} = savename;
end

end
