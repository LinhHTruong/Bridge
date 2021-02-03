function no_columns = pier_count_column(ptc_xyz, struct_threshold)
%% This function is to check the number of the column in data points of a column part of the pier
% Input:
%       ptc_xyz             : [Nx3] Point cloud of the column part
%       struct_threshold    : structural thresholds    
% Output:
%       no_columns          : the number of the columns
% Demo:
% ptc_xyz = pier_cluster_ptc_xyz;
%% Generate the elevation section
tic
ptc_min_z = min(ptc_xyz(:,3));
ptc_max_z = max(ptc_xyz(:,3));
num_sects = ceil((ptc_max_z - ptc_min_z)/struct_threshold.section_interval);
sects_locs = linspace(ptc_min_z, ptc_max_z, num_sects);

% Generate 2D cells
cluster_count = inf(numel(sects_locs)-1,1);
for i=1:numel(sects_locs)-1
    % Retrieve the data in the section
    mask = (sects_locs(i) <= ptc_xyz(:,3))&(ptc_xyz(:,3) <= sects_locs(i+1));
    sect_ptc_xyz = ptc_xyz(mask,1:3);
    
    % Generate 2D cell
    cell_tree = OctQuadtree(sect_ptc_xyz,'max_size', 0.5*struct_threshold.section_width);
    leaf_cell_ids = Node_Leaf(cell_tree);
    
    % Connectivity cluster
    region_info = cell_connectivity_segmentation(cell_tree, leaf_cell_ids);
    % Determine the number of the cluster
    cluster_count(i) = max(region_info(:,2));
    clear mask sect_ptc_xyz cell_tree leaf_cell_ids region_info
end

% Count coocur the number of clusters
count = unique(cluster_count);
region_count = histcounts(cluster_count, max(count));
region_count = region_count(count);
region_stat(:, 1) = unique(count);
region_stat(:, 2) = region_count;

% Determine the number of the column
[~, max_id] = max(region_stat(:, 2));
no_columns = region_stat(max_id,1);
