function [COM_new,I_new, m_new] = AdaptSegmentCOMAndInetia(COM_or,I_or,m_or,...
    COM_mAdded,I_mAdded, m_mAdded)
%AdaptSegmentCOMAndInetia Summary of this function goes here
%   Detailed explanation goes here


% new location COM
COM_new = (COM_or*m_or + COM_mAdded*m_mAdded)./(m_or + m_mAdded);

% Inertia x-axis
dCOM = sqrt((COM_or(2)-COM_new(2)).^2 +(COM_or(3)-COM_new(3)).^2);
dMass = sqrt((COM_mAdded(2)-COM_new(2)).^2 +(COM_mAdded(3)-COM_new(3)).^2);
I_new_x = I_or(1) + m_or * dCOM.^2 + I_mAdded(1) + m_mAdded * dMass.^2;


% Inertia y-axis
dCOM = sqrt((COM_or(1)-COM_new(1)).^2 +(COM_or(3)-COM_new(3)).^2);
dMass = sqrt((COM_mAdded(1)-COM_new(1)).^2 +(COM_mAdded(3)-COM_new(3)).^2);
I_new_y = I_or(2) + m_or * dCOM.^2 + I_mAdded(2) + m_mAdded * dMass.^2;


% Inertia z-axis
dCOM = sqrt((COM_or(1)-COM_new(1)).^2 +(COM_or(2)-COM_new(2)).^2);
dMass = sqrt((COM_mAdded(1)-COM_new(1)).^2 +(COM_mAdded(2)-COM_new(2)).^2);
I_new_z = I_or(3) + m_or * dCOM.^2 + I_mAdded(3) + m_mAdded * dMass.^2;

% new intertia
I_new = [I_new_x I_new_y I_new_z];

% new total mass
m_new = m_or + m_mAdded;



end