function [mtot] = getModelMass(ModelFile)
% --------------------------------------------------------------------------
% getModelMass
%   Computes to sum of all segment masses of an opensim model. 
%   Copied from https://github.com/KULeuvenNeuromechanics/NeuromechanicsToolkit/blob/75af1dead6960936d323b862b463ff789fb08692/OpenSimAPI/getModelMass.m
% 
% INPUT:
%   - Modelfile -
%   * path to the opensim model
%
% OUTPUT:
%   - mtot -
%   * total mass of the opensim model
% 
% Original author: Maarten Afschrift
% Original date: 6/Dec/2021
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import org.opensim.modeling.*
m = Model(ModelFile);
mtot = 0;

for i=1:m.getBodySet.getSize()
    mtot = mtot + m.getBodySet.get(i-1).getMass();
end




end