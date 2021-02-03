function [cells_plane, cells_nonplane, cells_region] = bridge_superstructure_surface_extraction(Tree, cells_plane, cells_nonplane, threshold, surface_position, running_time)
%% This function is to extract the top and bottom surfaces of the superstructure based on a cell base region growing
% Input:
% Output:
% Demo:
% Tree = OQTR;
% cells_plane = Cells_Plane
% cells_nonplane = Cells_nonPlane
% threshold = THRESHOLD
% surface_position = 'top'
% running_time = true

% Region growing
[region_info, cells_plane] = cell_region_growing_segmentation(Tree, cells_plane, threshold, surface_position, running_time);
% Back and Forward filtering: For the cells are shared by multiple segments/regions
cells_region = cell_backward_forward_filtering(Tree, cells_plane, region_info, threshold, running_time);
%
% Forward region growing: for the cells on the boundary of the region to
% searching points out off the peak (the surface) of the cells out off the region

[cells_plane, cells_region] = cell_forward_region_growing(Tree, cells_plane, cells_region, threshold, running_time);

% Back and Forward filtering: For the cells are shared by multiple segments/regions
[cells_plane, cells_region] = cell_boundary_filtering(Tree, cells_plane, cells_region, threshold, running_time);

% Forward searching the points on non-planar cells
[cells_nonplane, cells_region] = cell_non_plane_filtering(Tree, cells_nonplane, cells_region, threshold, running_time);

% Merging regions
cells_region = segment_merging(Tree, cells_region, threshold, running_time);

