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
% Last edit by:
% Last edit date:
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