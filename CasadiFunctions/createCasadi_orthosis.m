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


%% Define variables
% all inputs
qk = SX.sym('q',n_coord,1);
qdotk = SX.sym('qdot',n_coord,1);
GRFk = SX.sym('GRF',6);

% all outputs
Tk = SX(n_coord,1);

if ~isempty(S.orthosis)
    for i=1:length(S.orthosis.settings)
        otype = S.orthosis.settings{i}.type;
        switch otype
            case 'passive'
                f_passiveOrthosis_i = createCasadi_passiveOrthosisDescriptions(S,model_info,S.orthosis.settings{i});
                Tk_i = f_passiveOrthosis_i(qk,qdotk,GRFk);
                Tk = Tk + Tk_i;
    
            otherwise
                error(['Unknown orthosis type "' otype '".'])
        end
    end
end


%%

f_orthosis = Function('f_orthosis',{qk,qdotk,GRFk},{Tk},{'qk','qdotk','GRFk'},{'Tk'});



end