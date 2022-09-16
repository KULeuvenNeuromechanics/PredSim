function [f_orthosis] = createCasadi_orthosis(S,model_info)
% --------------------------------------------------------------------------
% createCasadi_orthosis
%   This function creates a casadifunction to calculate the total orthosis
%   torque on every joint.
% 
% INPUT:
%   - S -
%   * setting structure S
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
%
% OUTPUT:
%   - f_orthosis -
%   * casadi function to calculate orthosis torque
% 
% Original author: Lars D'Hondt
% Original date: 14/September/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.*

n_coord = model_info.ExtFunIO.jointi.nq.all;
N = S.solver.N_meshes;

%% Define variables
% all inputs
qk = SX.sym('q',n_coord,1);
qdotk = SX.sym('qdot',n_coord,1);

% all outputs
Tk = SX('T',n_coord,1);


for i=1:length(S.orthosis.settings)
    otype = S.orthosis.settings{i}.type;
    switch otype
        case 'passive'
            f_passiveOrthosis_i = createCasadi_passiveOrthosisDescriptions(S,model_info,S.orthosis.settings{i});
            Tk_i = f_passiveOrthosis_i(qk,qdotk);
            Tk = Tk + Tk_i;
    end
end





















