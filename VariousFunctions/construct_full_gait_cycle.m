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
            error("No symmetry defined for '%s'",field)
    end


    % 0-50%
    var1stHalfGC = res.(field);

    % 50-100%
    var2ndHalfGC = var1stHalfGC;
    var2ndHalfGC(idxInvA,:) = var2ndHalfGC(idxInvB,:);
    var2ndHalfGC(idxOpp,:) = -var2ndHalfGC(idxOpp,:);

    if strcmp(field,'q') % offset forward position
        var2ndHalfGC(model_info.ExtFunIO.jointi.base_forward,:) = ...
            var2ndHalfGC(model_info.ExtFunIO.jointi.base_forward,:) + ...
            var1stHalfGC(model_info.ExtFunIO.jointi.base_forward,end) - ...
            var1stHalfGC(model_info.ExtFunIO.jointi.base_forward,1);
    end

    res_GC.(field) = [var1stHalfGC, var2ndHalfGC];

end

end % end of function