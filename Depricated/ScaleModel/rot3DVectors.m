function rotated = rot3DVectors(rot, vecTrajs)
% Rotate any N number of 3D points/vectors
% USAGE: rotated = rot3DVectors(rot, vecTrajs)
%        rot is 3x3 rotation matrix
%        vecTrajs, Matrix of 3D trajectories (i.e. ntime x 3N cols)
% Ajay Seth

[nt, nc] = size(vecTrajs);

if rem(nc,3),
    error('Input trajectories must have 3 components each.');
end

for I = 1:nc/3,
    vecTrajs(:,3*I-2:3*I) = [rot*vecTrajs(:,3*I-2:3*I)']';
end

rotated = vecTrajs;