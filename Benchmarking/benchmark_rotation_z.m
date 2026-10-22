function rotation = benchmark_rotation_z(angle)
% Three-dimensional Z-axis rotation; angle is in radians.
rotation = [cos(angle), -sin(angle), 0; sin(angle), cos(angle), 0; 0, 0, 1];
end
