function scaleContactSpheres(osim_path_in,osim_path_out,scale)
% This code scales:
% 1. the radius and location of the contact spheres based on the scaling 
% factor of the foot. 
% 2. the stiffness of the contact sphere
% 3. the dissipation of the contact sphere
% 
% INPUT:
%
%   - osim_path_in -
%   * Path to the model whose contact spheres you want to scale
%
%   - osim_path_out -
%   * Path to where the model with scaled contact spheres will be written.
% 
%   - scale -
%   * a struct with the fields 'stiffness', 'dissipation', 'foot_left',
%       and 'foot_right' containing the scaling factors to apply.
% 
% OUTPUT:
%   - a new model with the scaled contact spheres
% 
% Original author: Bram Van Den Bosch
% Original date: 27/March/2023
% --------------------------------------------------------------------------

%% 1. get contact spheres

contact_spheres = getContactSpheres(osim_path_in);

%Get the position of toes origin with respect to calcaneus
import org.opensim.modeling.*;


%% 2. calculate new values

for i=1:length(contact_spheres)
    contact_spheres(i).stiffness = contact_spheres(i).stiffness(:,:)*scale.stiffness;
    contact_spheres(i).dissipation = contact_spheres(i).dissipation*scale.dissipation;

    name = contact_spheres(i).name;
    [leftname,~] = mirrorName(name);
    if strcmp(name,leftname)
        contact_spheres(i).radius = contact_spheres(i).radius*scale.foot_left(1);
        contact_spheres(i).location = contact_spheres(i).location.*scale.foot_left;
    else
        contact_spheres(i).radius = contact_spheres(i).radius*scale.foot_right(1);
        contact_spheres(i).location = contact_spheres(i).location.*scale.foot_right;
    end
end

%% 3. write values to model

model = Model(osim_path_in);
model.finalizeConnections();

forceSet = model.getForceSet();
geometrySet = model.getContactGeometrySet();

for i=0:forceSet.getSize()-1
    force_i = forceSet.get(i);
    if strcmp(force_i.getConcreteClassName(),"SmoothSphereHalfSpaceForce")

        force_i = SmoothSphereHalfSpaceForce.safeDownCast(force_i);

        socket1 = force_i.getSocket("sphere");
        socket1_objPath = char(socket1.getConnecteePath());
        [~,socket1_objName,~] = fileparts(socket1_objPath);
        idx = find(strcmp(socket1_objName,{contact_spheres(:).name}));
        geo = geometrySet.get(socket1_objName);
        geo1 = ContactSphere.safeDownCast(geo);

        % set stiffness
        force_i.set_stiffness(contact_spheres(idx).stiffness)
        % set dissipation
        force_i.set_dissipation(contact_spheres(idx).dissipation)
        % set radius
        geo1.setRadius(contact_spheres(idx).radius);
        % set location
        contactSphereLoc = Vec3.createFromMat(contact_spheres(idx).location);
        geo1.setLocation(contactSphereLoc);

    end
end

model.finalizeConnections();
model.initSystem();
model.print(osim_path_out);

end