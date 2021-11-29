% This scripts implements the Hunt-Crossley contact model described in
% Sherman et al. (2011). We use a smooth approximation in order to make the 
% model twice continuously differentiable.
%
% INPUTS:
%   - stiffness of the contact sphere
%   - radius of the contact sphere
%   - dissipation: coefficient of effective dissipation 
%   - normal: orientiation of the vector perpendicular to the contact
%       sphere; [0,1,0] for flat ground and perpendicular contact
%   - transitionVelocity: velocity at which static friction reaches its 
%       peak value
%   - staticFriction: coefficient of static friction 
%   - dynamicFriction: coefficient of dynamic friction 
%   - viscousFriction: coefficient of viscous friction 
%   - spherePos_inB: position of the contact sphere origin in the body
%       frame
%   - posFrame_inG: position of the body frame origin in ground 
%        (OpenSim: getPositionInGround)
%   - linVelFrame_inG: linear velocity of the body frame in ground
%       (OpenSim: getVelocityInGround [1])
%   - angVelFrame_InG: angular velocity of the body frame in ground
%       (OpenSim: getVelocityInGround [0])
%   - rotBtoG_inG: rotational component of the transform from the body  
%       frame to the ground frame in ground (OpenSim: getBodyTransform.R)
%   - trBtoG_inG: translational component of the transform from the body  
%       frame to the ground frame in ground (OpenSim: getBodyTransform.T)
%
% OUTPUT:
%   - force: applied contact force
%
% Author: Antoine Falisse
% Date: 1/7/2019
%
function force = HCContactModel(stiffness,radius,dissipation,normal,...
    transitionVelocity,staticFriction,dynamicFriction,viscousFriction,...
    spherePos_inB,posFrame_inG,linVelFrame_inG,angVelFrame_InG,...
    rotBtoG_inG,trBtoG_inG) 

% Reshape rotation matrix
Rot_l = (reshape(rotBtoG_inG,3,3))';
% Express sphere position in ground
spherePos = (Rot_l*spherePos_inB+trBtoG_inG)';
% Contact point position
temp = spherePos-[0,radius,0];
indentation = -temp(2);
indentation_vec = -temp(2)*normal;
% Contact point velocity
v = linVelFrame_inG' + cross(angVelFrame_InG',spherePos - [0,radius,0] +...
    0.5*indentation_vec - posFrame_inG');
vnormal = v(1)*normal(1) + v(2)*normal(2) + v(3)*normal(3);
vtangent = v - vnormal*normal;
indentationVel = -vnormal;
% Constant values
eps = 1e-5;
eps2 = 1e-16;
bv = 50;
bd = 300; 
% Stiffness force
k = 0.5*(stiffness)^(2/3);
fH = (4/3)*k*sqrt(radius*k)*((sqrt(indentation*indentation+eps))^(3/2));
% Dissipation force
fHd = fH*(1+1.5*dissipation*indentationVel); 
fn = (0.5*tanh(bv*(indentationVel+1/(1.5*dissipation)))+0.5 + eps2)*...
    (0.5*tanh(bd*indentation)+0.5 + eps2)*fHd;
force = fn.*normal;
% Friction force
aux = (vtangent(1)).^2 + (vtangent(2)).^2 + (vtangent(3)).^2 + eps; 
vslip = aux.^(0.5);
vrel = vslip/transitionVelocity;
ffriction = fn*(min(vrel,1)*(dynamicFriction + 2 * (staticFriction - ...
    dynamicFriction) / (1 + vrel*vrel)) + viscousFriction*vslip);
% Contact force
force = force + ffriction*(-vtangent) / vslip;

end

