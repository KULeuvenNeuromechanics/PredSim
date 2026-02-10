function scaleContactSpheres(genericModelPath,modelpath_in,modelout_path,scale)
% This code scales:
% 1. the radius of the contact spheres based on the scaling factor of the
% foot.
% 2. the stiffness of the contact sphere
% 3. the dissipation of the contact sphere
%
% INPUT:
%   - genericModelPath: Path to the generic model.
%
%   - modelpath_in: Path to the model whose contact spheres you want to scale
%
%   - modelout_path: Path to where the model with scaled contact spheres
%   will be written
%
%   - scale: a struct with optional fields:
%       * 'stiffness': scaling factor for stiffness (default: 1.0)
%       * 'dissipation': scaling factor for dissipation (default: 1.0)
%       * 'factorR': scaling factor for right foot radius (default: auto-calculated)
%       * 'factorL': scaling factor for left foot radius (default: auto-calculated)
%
% OUTPUT:
%   - a new model with the scaled contact spheres
%
% Original author: Bram Van Den Bosch
% Original date: 27/March/2023
%
% Last edit by: Ellis Van Can
% Last edit date: 13/November/2025
% --------------------------------------------------------------------------

%% Handle optional inputs and set defaults
% Default scale struct if not provided or empty
if isempty(scale)
    scale = struct();
end

% Set default values for missing fields
if ~isfield(scale, 'stiffness')
    scale.stiffness = 1.0;
end
if ~isfield(scale, 'dissipation')
    scale.dissipation = 1.0;
end

%% Get contact spheres
contact_spheres = get_contact_spheres(modelpath_in);

%% Get the position of toes origin with respect to calcaneus
import org.opensim.modeling.*;
model_1 = Model(genericModelPath);
toesInCalcn_r_1 = model_1.getJointSet.get('mtp_r').get_frames(0).get_translation.getAsMat';
toesInCalcn_l_1 = model_1.getJointSet.get('mtp_l').get_frames(0).get_translation.getAsMat';

model_2 = Model(modelpath_in);
toesInCalcn_r_2 = model_2.getJointSet.get('mtp_r').get_frames(0).get_translation.getAsMat';
toesInCalcn_l_2 = model_2.getJointSet.get('mtp_l').get_frames(0).get_translation.getAsMat';

% Determine scaling factors for radius
if isfield(scale, 'factorR') && ~isempty(scale.factorR)
    factorR = scale.factorR;
else
    factorR = sqrt(sum(toesInCalcn_r_2.^2))/sqrt(sum(toesInCalcn_r_1.^2));
end

if isfield(scale, 'factorL') && ~isempty(scale.factorL)
    factorL = scale.factorL;
else
    factorL = sqrt(sum(toesInCalcn_l_2.^2))/sqrt(sum(toesInCalcn_l_1.^2));
end

%% Calculate new values
for i = 1:length(contact_spheres)
    contact_spheres(i).stiffness = contact_spheres(i).stiffness * scale.stiffness;
    contact_spheres(i).dissipation = contact_spheres(i).dissipation * scale.dissipation;
    
    if contains(contact_spheres(i).body, '_r')
        contact_spheres(i).radius = contact_spheres(i).radius * factorR;
    else
        contact_spheres(i).radius = contact_spheres(i).radius * factorL;
    end
end

%% Write values to model
model = Model(modelpath_in);
model.finalizeConnections();
forceSet = model.getForceSet();
geometrySet = model.getContactGeometrySet();

for i = 0:forceSet.getSize()-1
    force_i = forceSet.get(i);
    
    if strcmp(char(force_i.getConcreteClassName()), 'SmoothSphereHalfSpaceForce')
        force_i = SmoothSphereHalfSpaceForce.safeDownCast(force_i);
        socket1 = force_i.getSocket('sphere');
        socket1_objPath = char(socket1.getConnecteePath());
        [~, socket1_objName, ~] = fileparts(socket1_objPath);
        
        idx = find(strcmp(socket1_objName, {contact_spheres(:).name}));
        
        if ~isempty(idx)
            geo = geometrySet.get(socket1_objName);
            geo1 = ContactSphere.safeDownCast(geo);
            
            % Set stiffness
            force_i.set_stiffness(contact_spheres(idx).stiffness);
            
            % Set dissipation
            force_i.set_dissipation(contact_spheres(idx).dissipation);
            
            % Set radius
            geo1.setRadius(contact_spheres(idx).radius);
        end
    end
end

model.finalizeConnections();
model.initSystem();
model.print(modelout_path);
end