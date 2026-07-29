function [bounds] = scaleBounds(S,model_info,bounds_nsc,scaling)
% --------------------------------------------------------------------------
% scaleBounds
%   This script scales bounds for the scaled optimisation variables.
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
%   - bounds_nsc -
%   * boundaries for all optimisation variables (not scaled)
%
%   - scaling -
%   * scale factors for all optimisation variables
%
% OUTPUT:
%   - bounds -
%   * boundaries for all scaled optimisation variables
%
% 
% Original author: Lars D'Hondt
% Original date: 6/April/2023
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

bounds = bounds_nsc;

% Qs
bounds.Qs.lower = (bounds_nsc.Qs.lower)./scaling.Qs;
bounds.Qs.upper = (bounds_nsc.Qs.upper)./scaling.Qs;
% Qdots
bounds.Qdots.lower = (bounds_nsc.Qdots.lower)./scaling.Qdots;
bounds.Qdots.upper = (bounds_nsc.Qdots.upper)./scaling.Qdots;
% Qdotdots
bounds.Qdotdots.lower = (bounds_nsc.Qdotdots.lower)./scaling.Qdotdots;
bounds.Qdotdots.upper = (bounds_nsc.Qdotdots.upper)./scaling.Qdotdots;
bounds.Qdotdots.lower(isnan(bounds.Qdotdots.lower)) = 0; % not sure when this check is needed, so I'm keeping it
bounds.Qdotdots.upper(isnan(bounds.Qdotdots.upper)) = 0;
% Muscle-tendon forces
bounds_nsc.FTtilde.lower    = (bounds_nsc.FTtilde.lower)./scaling.FTtilde;
bounds_nsc.FTtilde.upper    = (bounds_nsc.FTtilde.upper)./scaling.FTtilde;

% We impose the initial position of pelvis_tx to be 0
bounds.Qs_0.lower = bounds.Qs.lower;
bounds.Qs_0.upper = bounds.Qs.upper;
bounds.Qs_0.lower(model_info.ExtFunIO.jointi.base_forward) = 0;
bounds.Qs_0.upper(model_info.ExtFunIO.jointi.base_forward) = 0;
bounds.Qs_0.lower(model_info.ExtFunIO.jointi.base_lateral) = 0;
bounds.Qs_0.upper(model_info.ExtFunIO.jointi.base_lateral) = 0;


end % end of function