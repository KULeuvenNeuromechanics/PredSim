function [S, osimConstraints] = kinematic_constraints_helper(S, osim_path)
% --------------------------------------------------------------------------
% kinematic_constraints_helper
%   Detect kinematic constraints in the opensim model and update
%   S.OpenSimADOptions with the required inputs and outputs.
% 
%
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%
% OUTPUT:
%   - S -
%   * setting structure S
% 
% Original author: Lars D'Hondt
% Original date: 24 September 2025
% --------------------------------------------------------------------------

import org.opensim.modeling.*;

model = Model(osim_path);
model.initSystem;
constr_set = model.getConstraintSet();

osimConstraints = [];

for i=0:constr_set.getSize()-1

    constr_i = constr_set.get(i);
    constr_name_i = char(constr_i.getName());

    % conditionally enforce constraints?
%     isEnforced_i = constr_i.get_isEnforced(); 
    
    if strcmp(constr_i.getConcreteClassName,'PointConstraint')
    
        pconstr = PointConstraint.safeDownCast(constr_i);
        
        con_i_body(1).body = char(pconstr.getConnectee('body_1').getName());
        con_i_body(2).body = char(pconstr.getConnectee('body_2').getName());
    
        con_i_body(1).point_in_body = pconstr.get_location_body_1().getAsMat();
        con_i_body(2).point_in_body = pconstr.get_location_body_2().getAsMat();

        con_i_body(1).name = ['osimConstraint_',constr_name_i,'_1'];
        con_i_body(2).name = ['osimConstraint_',constr_name_i,'_2'];
        
        con_i_body(1).reference_frame = 'ground';
        con_i_body(2).reference_frame = 'ground';
            % use last common 'parent' instead? does this matter a lot?

        S.OpenSimADOptions.export3DPositions = ...
            [S.OpenSimADOptions.export3DPositions(:); con_i_body(:)];

        S.OpenSimADOptions.input3DBodyForces = ...
            [S.OpenSimADOptions.input3DBodyForces(:); con_i_body(:)];
    
        osimConstraints{end+1} = constr_name_i;
    end
end




end % end of function
