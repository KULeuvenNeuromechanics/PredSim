function [gomenuka2014] = adapt_model_Gomenuka(S,osim_path)
%adapt_model_Gomenuka Adapts default .osim to replicate experiments of
%gomenuka 2014
%   input arguments:
%       (1) S: default settings structure
%       (2) S: osim_path
%   output arguments:
%       (1) output structure with path information about the exported
%       models
import org.opensim.modeling.*;
slopes = [0 0.07 0.15 0 0.07 0.15];
MassTorso = [zeros(1,3), zeros(1,3)+0.25];
ModelOut = {'_default'; '_slope07'; '_slope15';...
    '_torso125pct'; '_slope07_torso125pct'; '_slope15_torso125pct'};

model_mass = getModelMass(osim_path);
ct = 1;
for imodels = 1:length(ModelOut)

    % model prop
    slope = slopes(imodels);
    added_mass = model_mass * MassTorso(imodels);

    % adapt gravity vector of opensim model
    % read model
    modSel = Model(osim_path);    

    % adapt gravity based on slope
    g = modSel.getGravity();
    gv = [g.get(0) g.get(1) g.get(2)];
    fi = atan(slope);
    R = rotz(fi);
    R = R(1:3, 1:3);
    grav_tilt = gv*R;
    gravSlope = Vec3(grav_tilt(1), grav_tilt(2), grav_tilt(3));
    modSel.setGravity(gravSlope);

    % adapt torso mass
    if added_mass ~= 0 
        BodySel = modSel.getBodySet().get('torso');
        mBodySel = BodySel.getMass();
        BodySel.setMass(mBodySel + added_mass);
    end

        % init the model (not sure if this is needed)
    modSel.initSystem();
    
    % save model
    model_name = [S.subject.name ModelOut{imodels}];
    out_folder =  fullfile(S.misc.main_path,'Subjects',model_name);
    out_modelname = fullfile(out_folder,[model_name '.osim']);
    if ~isfolder(out_folder)
        mkdir(out_folder)
    end
    modSel.print(out_modelname);

    % update output structure
    gomenuka2014.modelnames{ct} = model_name;
    gomenuka2014.osim_path{ct} = out_modelname;
    ct = ct+1;

    % save 
    clear mSel;

end

end