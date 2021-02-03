function pier_cluster = bridge_pier_type(comp_ptc_xyz, struct_threshold)
%% This function is to classify a type of the pier
% Input:
%       comp_ptc_xyz            : [Nx3] - x, y, z coordinates of the components
%       struct_threshold        :
% Output:
%       pier_cluster            : [Nx2] - N is the number of the cluster, [s_elev, e_elev]
% 
% Demo:

%% Generate the cross-section
ptc_min_z = min(comp_ptc_xyz(:,3));
ptc_max_z = max(comp_ptc_xyz(:,3));
num_sects = ceil((ptc_max_z - ptc_min_z)/struct_threshold.section_interval);
sects_locs = linspace(ptc_min_z, ptc_max_z, num_sects);

%% Compute the cross-section dimension
sects_dim = inf(numel(sects_locs)-1, 4);%[start_z, end_z, short_edge, long_edge]
for i = 1:size(sects_dim)
    % Retrieve the points in the section
    mask = (sects_locs(i) <= comp_ptc_xyz(:,3))&(comp_ptc_xyz(:,3) <= sects_locs(i+1));
    sect_ptc_xyz = comp_ptc_xyz(mask,1:3);
    
    % Compute the dimensions throung minimum bounding box
    rot_sect_mbb = min2DBoundingBox(sect_ptc_xyz(:,1:2)');
    sects_dim(i,:) = [sects_locs(i), sects_locs(i+1), rot_sect_mbb.short_edge, rot_sect_mbb.long_edge];
end

%% Clustering the sects based on its length
cluster_ids = dbscan(sects_dim(:,4), 0.5*struct_threshold.section_length,4);
num_cluster = max(cluster_ids);
% Preallocation for a cluster 
pier_cluster = [];
for i = 1:num_cluster
    % Retrieve the section in the cluster
    mask = cluster_ids == i;
    cluster_z = sects_dim(mask,:);
    
    % Use BD scan to cluster the group according to elevation
    cluster_z_ids = dbscan(cluster_z(:,1), 3.0*struct_threshold.section_interval,2);
    
    for j = 1:max(cluster_z_ids)
        % Retrieve locations of each cluster
        mask = cluster_z_ids == j;
        pier_cluster = [pier_cluster; [min(cluster_z(mask,1)), max(cluster_z(mask,2)), mean(cluster_z(mask,[3,4]))]];
    end
end
%% Classified the cluster
pier_component_code = inf(size(pier_cluster,1),1);
[~, mask] = sort(pier_cluster(:,1), 'descend');
pier_cluster = pier_cluster(mask,:);
% Find the id of the column
pier_cluster_height = pier_cluster(:,2) - pier_cluster(:,1);
[~, max_id] = max(pier_cluster_height);

% Determine the pier cap
mask = pier_cluster(:,2) > mean(pier_cluster(max_id,[1,2]));
pier_component_code(mask) = 1; % Pier cap

% Determine the pile cap
mask = pier_cluster(:,1) < mean(pier_cluster(max_id,[1,2]));
pier_component_code(mask) = 2; % Pile cap

% Assign the column
pier_component_code(max_id) = 3; %3 coloum 

% Assembly
pier_cluster = [pier_cluster, pier_component_code];

    