function [ptc_local_ids, ptc_xyz] = extract_ptc_2_planes(ptc_xyz, st_plane, nd_plane)
%% This fuction is to extract the points between two planes
% Input:
%   ptc_xyz             : [Nx6] x, y, z
%   st_plane            : [1x6] x, y, z, nx, ny, nz
%   st_plane            : [1x6] x, y, z, nx, ny, nz

% Output:
%   ptc_local_ids       : [Kx1] a local indices of the point between the planes
%   ptc_xyz             : [Kx3] x, y, z
%
% Demo:
% ptc_xyz = hor_seg_ptc_xyz;
% st_plane
% nd_plane

%% Adjust a direction of the normal
if dot(st_plane(4:6), nd_plane(4:6)) < 0
    nd_plane(4:6) = - nd_plane(4:6);
end

%% Compute the distance
dist_ptc_st_plane = dist_3Dpoints_3Dplane(ptc_xyz, st_plane);
dist_ptc_nd_plane = dist_3Dpoints_3Dplane(ptc_xyz, nd_plane);

%% Filter the points between two plane
mask = sign(dist_ptc_st_plane).*sign(dist_ptc_nd_plane) < 0;
ptc_local_ids = find(mask == true);
ptc_xyz = ptc_xyz(mask,1:3);