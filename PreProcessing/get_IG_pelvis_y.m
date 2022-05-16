function [model_info] = get_IG_pelvis_y(S,osim_path,model_info)
% --------------------------------------------------------------------------
% get_IG_pelvis_y
%   This function returns the pelvis height for quasi-random initial guess
%
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%       
% 
% Original author: Lars D'Hondt
% Original date: 12/April/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

if isempty(S.subject.IG_pelvis_y)
    import org.opensim.modeling.*;
    
    model = Model(osim_path);
    state = model.initSystem;
    
    calcn_or = model.getBodySet().get('calcn_r').findBaseFrame().getPositionInGround(state).getAsMat;
    
    model_info.IG_pelvis_y = -calcn_or(2)+25e-2;

else
    model_info.IG_pelvis_y = S.subject.IG_pelvis_y;

end









