function pt = cal_line_plane_intersection(line, plane)
%% This function is to determine an intersection points between the 3D
% line and 3D planes
% Input
% Output
% Ref:
%   http://geomalgorithms.com/a05-_intersect-1.html    
% Demo
%% Retrieve components
line_pt = line(1:3);
line_vect = line(4:6);
plane_pt = plane(1:3);
plane_vect = plane(4:6);

%% Calculate the intersection
s = dot(plane_vect, (plane_pt - line_pt))/dot(plane_vect, line_vect);

pt = line_pt + s*line_vect;


