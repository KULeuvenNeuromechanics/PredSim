function [] = removeContactSpheres(osim_path)
% --------------------------------------------------------------------------
% removeContactSpheres
%   This functions removes the contact spheres (SmoothSphereHalfSpaceForce 
%   and associated contact geometries) from a given .osim model.
% 
% INPUT:
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%
% Original author: Lars D'Hondt
% Original date: 3 January 2025
% --------------------------------------------------------------------------


%% load model
import org.opensim.modeling.*;
model = Model(osim_path);
model.finalizeConnections();

forceSet = model.getForceSet();
geometrySet = model.getContactGeometrySet();


%% remove contact spheres
contact_geometry = {};

% loop over forces to find contact forces
for i=0:forceSet.getSize()-1
    force_i = forceSet.get(i);
    if strcmp(force_i.getConcreteClassName(),"SmoothSphereHalfSpaceForce")

        force_i = SmoothSphereHalfSpaceForce.safeDownCast(force_i);

        % add name of half space to list of geometries to be removed
        socket1 = force_i.getSocket("half_space");
        socket1_objPath = char(socket1.getConnecteePath());
        [~,socket1_objName,~] = fileparts(socket1_objPath);
        contact_geometry{end+1} = socket1_objName;

        % add name of sphere to list of geometries to be removed
        socket1 = force_i.getSocket("sphere");
        socket1_objPath = char(socket1.getConnecteePath());
        [~,socket1_objName,~] = fileparts(socket1_objPath);
        contact_geometry{end+1} = socket1_objName;
  
        % remove force
        force_i.delete();
    end
end

% in case there were no contact spheres hence no geometries to remove
if isempty(contact_geometry)
    return
end

% remove contact geometries
for s=string(unique(contact_geometry))
    geometrySet.get(s).delete();

end



%% save model
model.finalizeConnections();
model.initSystem();
model.print(osim_path);


end % end of function
