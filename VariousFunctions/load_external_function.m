function [F] = load_external_function(S)
% --------------------------------------------------------------------------
% load_external_function
%   The external function performs inverse dynamics through the
%   OpenSim/Simbody C++ API. This external function is compiled as a dll from
%   which we create a Function instance using CasADi in MATLAB.
% 
% INPUT:
%   - S -
%   * setting structure S
%
% OUTPUT:
%   - F -
%   * function handle to the external function
% 
% Original author: Lars D'Hondt
% Original date: 10/May/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.*
pathmain = pwd;
cd(S.misc.subject_path)
F  = external('F',S.misc.external_function);
cd(pathmain);

