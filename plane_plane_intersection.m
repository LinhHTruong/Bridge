function [intersection_line, check] = plane_plane_intersection(plane_A, plane_B, varargin)
%plane_intersect computes the intersection of two planes(if any)
% Inputs: 
%       plane_A: [x, y, z, nx, ny, nz] center and normal of Plane 1
%       plane_B: [x, y, z, nx, ny, nz] center and normal of Plane 2
%
% Outputs:
%   P    is a point that lies on the interection straight line.
%   N    is the direction vector of the straight line
% check is an integer (0:Plane 1 and Plane 2 are parallel' 
%                              1:Plane 1 and Plane 2 coincide
%                              2:Plane 1 and Plane 2 intersect)
%
% Example:
% Determine the intersection of these two planes:
% 2x - 5y + 3z = 12 and 3x + 4y - 3z = 6
% The first plane is represented by the normal vector N1=[2 -5 3]
% and any arbitrary point that lies on the plane, ex: A1=[0 0 4]
% The second plane is represented by the normal vector N2=[3 4 -3]
% and any arbitrary point that lies on the plane, ex: A2=[0 0 -2]
%[P,N,check]=plane_intersect([2 -5 3],[0 0 4],[3 4 -3],[0 0 -2]);
%This function is written by :
%                             Nassim Khaled
%                             Wayne State University
%                             Research Assistant and Phd candidate

if numel(varargin) == 1
    max_normal_angle = varargin{:};
else
    max_normal_angle = 2.0;%10^-7;
end

% Preallocation a point on an intersection line
P = [0 0 0];

% Retrieve para.s of the plane
plane_A_pt = plane_A(1:3);
plane_A_norm = plane_A(4:6);

plane_B_pt = plane_B(1:3);
plane_B_norm = plane_B(4:6);

% Compute angle between two normals
cross_A_B = cross(plane_A_norm,plane_B_norm);
%  test if the two planes are parallel
if abs(dot(plane_A_norm, plane_B_norm)) >= cos(deg2rad(1.0*max_normal_angle))          % Plane 1 and Plane 2 are near parallel norm(cross) >= xxx
    v = plane_A_pt - plane_B_pt;
    if (dot(plane_A_norm,v) <= cos(deg2rad(90-max_normal_angle)))         
        check = -1;                    % Plane 1 and Plane 2 coincide
        intersection_line = [];
        return
    else 
        check = 0;                    % Plane 1 and Plane 2 are parallel
        intersection_line = [];
        return
    end
end

% The planes are intersected
check = 1;
 
% Plane 1 and Plane 2 intersect in a line first determine max abs coordinate of cross product
% maxc=find(abs(cross_A_B)==max(abs(cross_A_B)));
[~,maxc] = max(abs(cross_A_B));
%next, to get a point on the intersection line and zero the max coord, and solve for the other two
      
d1 = -dot(plane_A_norm, plane_A_pt);    % the constants in the Plane 1 equations
d2 = -dot(plane_B_norm, plane_B_pt);    % the constants in the Plane 2 equations

switch maxc
    case 1                   % intersect with x=0
        P(1)= 0;
        P(2) = (d2*plane_A_norm(3) - d1*plane_B_norm(3))/ cross_A_B(1);
        P(3) = (d1*plane_B_norm(2) - d2*plane_A_norm(2))/ cross_A_B(1);
    case 2                    % intersect with y=0
        P(1) = (d1*plane_B_norm(3) - d2*plane_A_norm(3))/ cross_A_B(2);
        P(2) = 0;
        P(3) = (d2*plane_A_norm(1) - d1*plane_B_norm(1))/ cross_A_B(2);
    case 3                    % intersect with z=0
        P(1) = (d2*plane_A_norm(2) - d1*plane_B_norm(2))/ cross_A_B(3);
        P(2) = (d1*plane_B_norm(1) - d2*plane_A_norm(1))/ cross_A_B(3);
        P(3) = 0;
end

% An intersection line
intersection_line = [P, cross_A_B/norm(cross_A_B)];
