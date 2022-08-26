function [res] = evaluate_external_function(S,model_info,Qs,Qdots,Qdotdots)
% --------------------------------------------------------------------------
% evaluate_external_function
%   Helper function to load and evaluate the external function. The external 
%   function performs inverse dynamics through the OpenSim/Simbody C++ API. 
%   This external function is compiled as a dll from which we create a 
%   Function instance using CasADi in MATLAB.
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - Qs -
%   * vector of coordinate angles/positions (rad, m)
%
%   - Qdots -
%   * vector of coordinate velocities (rad/s, m/s)
%
%   - Qdotdots -
%   * vector of coordinate accelerations (rad/s^2, m/s^2)
%
%   Alternative input: pass vector QsQdots (alternating Q and Qdot) instead
%   of Qs or Qdots, and give an empty input for Qdots or Qs respectively.
%
% OUTPUT:
%   - res -
%   * output vector of the external function. Variable type is MX or DM.
%
%
% Original author: Lars D'Hondt
% Original date: 28/Aug/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.*

% load external function
pathExt = fullfile(S.misc.subject_path,S.misc.external_function);
F  = external('F',pathExt);

% Adapt shape of input vectors
if isempty(Qs)
    QsQdots = vertcat(Qdots(:));
elseif isempty(Qdots)
    QsQdots = vertcat(Qs(:));
else
    QsQdots = MX(model_info.ExtFunIO.jointi.nq.all*2,1);
    QsQdots(1:2:end) = Qs(:);
    QsQdots(2:2:end) = Qdots(:);
end
Qdotdots = vertcat(Qdotdots(:));

% Set selected coordinate positions to zero to evaluate the external
% function (e.g. position of floating base in plane parallel to the ground)
if ~isempty(S.misc.coordPosZeroForExternal)
    for i=1:length(S.misc.coordPosZeroForExternal)
        if isfield(model_info.ExtFunIO.coordi,S.misc.coordPosZeroForExternal{i})
            idxToZero = model_info.ExtFunIO.coordi.(S.misc.coordPosZeroForExternal{i});
            QsQdots(idxToZero*2-1) = 0;
        end
    end
end

% Evaluate F
res = F([QsQdots; Qdotdots]);

end

