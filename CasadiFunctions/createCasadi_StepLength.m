function [f_getCalcnOriginInWorldFrame,f_getStepLength] = createCasadi_StepLength(S,model_info)
% --------------------------------------------------------------------------
% createCasadi_StepLength
%   Create CasADi functions to describe the position of the calcaneus
%   origin and the steplength.
%
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% OUTPUT:
%   - f_getCalcnOriginInWorldFrame -
%   * Origin of calcn_r and calcn_l bodies in world reference frame for the
%   given coordinate positions
%
%   - f_getStepLength -
%   * Step length for the given coordinate positions
%
% Original author: Lars D'Hondt
% Original date: 16/May/2022
%
% Last edit by:
% Last edit date:
% --------------------------------------------------------------------------

import casadi.*

qs = MX.sym('qs',model_info.ExtFunIO.jointi.nq.all);

% Create zero (sparse) input vector for external function
F_ext_input = MX(model_info.ExtFunIO.input.nInputs,1);
% Assign Qs
F_ext_input(model_info.ExtFunIO.input.Qs.all,1) = qs;

F = external('F',replace(fullfile(S.misc.subject_path,S.misc.external_function),'\','/'));

F_out = F(F_ext_input);

calcn_or_r = F_out(model_info.ExtFunIO.position.calcn_r);
calcn_or_l = F_out(model_info.ExtFunIO.position.calcn_l);


f_getCalcnOriginInWorldFrame = Function('f_getCalcnOriginInWorldFrame',{qs},...
    {calcn_or_r,calcn_or_l},{'q'},{'calcn_or_r','calcn_or_l'});

%%
qs_1 = MX.sym('qs_1',model_info.ExtFunIO.jointi.nq.all);
qs_end = MX.sym('qs_end',model_info.ExtFunIO.jointi.nq.all);

[calcn_or_r_1,calcn_or_l_1] = f_getCalcnOriginInWorldFrame(qs_1);
[calcn_or_r_end,calcn_or_l_end] = f_getCalcnOriginInWorldFrame(qs_end);

step_length_l = calcn_or_l_end(1) - calcn_or_l_1(1);
step_length_r = calcn_or_r_end(1) - calcn_or_r_1(1);

f_getStepLength = Function('f_getStepLength',{qs_1,qs_end},{step_length_r,step_length_l},...
    {'qs_1','qs_end'},{'step_length_r','step_length_l'});


end % end of function