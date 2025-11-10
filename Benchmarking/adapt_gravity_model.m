function [] = adapt_gravity_model(modelpath,slope,outpathmodel)
% adapt_gravity_model adapt gravity in opensim model


import org.opensim.modeling.*;
% read model
modSel = Model(modelpath);

% adapt gravity based on slope
g = modSel.getGravity();
gv = [g.get(0) g.get(1) g.get(2)];
fi = atan(slope/100);
R = rotz(fi);
R = R(1:3, 1:3);
grav_tilt = gv*R;
gravSlope = Vec3(grav_tilt(1), grav_tilt(2), grav_tilt(3));
modSel.setGravity(gravSlope);

% create output folder if needed
[outfolder,~,~]=fileparts(outpathmodel);
if ~isfolder(outfolder)
    mkdir(outfolder);
end

% save model
modSel.print(outpathmodel);



end