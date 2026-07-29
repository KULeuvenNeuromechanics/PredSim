function [res_GC] = construct_full_gait_cycle(model_info, res)
% --------------------------------------------------------------------------
% construct_full_gait_cycle
%   Constructs full gait cycle from half gait cycle based on symmetry
%   information.
%
% INPUT:
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
%   - res -
%   * Struct with variables for half gait cycle. Time is horizontal.
%
%
% OUTPUT:
%   - res_GC -
%   * Struct with variables for half gait cycle. Time is horizontal.
% 
% Original author: Lars D'Hondt
% Original date: 23 January 2024
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



for field = string(fieldnames(res)')

    % get symmetry indices
    switch field
        case 'q' 
            idxInvA = model_info.ExtFunIO.symQs.QsInvA;
            idxInvB = model_info.ExtFunIO.symQs.QsInvB;
            idxOpp = model_info.ExtFunIO.symQs.QsOpp;

        case {'qdot', 'qddot', 'T', 'M'}
            idxInvA = model_info.ExtFunIO.symQs.QdotsInvA;
            idxInvB = model_info.ExtFunIO.symQs.QdotsInvB;
            idxOpp = model_info.ExtFunIO.symQs.QsOpp;

        case {'a', 'vA', 'FTtilde', 'dFTtilde'}
            idxInvA = model_info.ExtFunIO.symQs.MusInvA;
            idxInvB = model_info.ExtFunIO.symQs.MusInvB;
            idxOpp = [];

        case {'a_a', 'e_a'}
            idxInvA = model_info.ExtFunIO.symQs.ActInvA;
            idxInvB = model_info.ExtFunIO.symQs.ActInvB;
            idxOpp = model_info.ExtFunIO.symQs.ActOpp;

        otherwise
            idxInvA = [];
            idxInvB = [];
            idxOpp = [];

    end


    % 0-50%
    var1stHalfGC = res.(field);

    % 50-100%
    var2ndHalfGC = nan(size(var1stHalfGC));
    var2ndHalfGC(idxInvA,:) = var1stHalfGC(idxInvB,:);
    var2ndHalfGC(idxOpp,:) = -var1stHalfGC(idxOpp,:);

    if strcmp(field,'q') % offset forward position
        var2ndHalfGC(model_info.ExtFunIO.jointi.base_forward,:) = ...
            var2ndHalfGC(model_info.ExtFunIO.jointi.base_forward,:) + ...
            var1stHalfGC(model_info.ExtFunIO.jointi.base_forward,end) - ...
            var1stHalfGC(model_info.ExtFunIO.jointi.base_forward,1);
    end

    res_GC.(field) = [var1stHalfGC, var2ndHalfGC];

end

end % end of function