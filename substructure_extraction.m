function subTree = substructure_extraction (Tree, active_pts, road_surface, struct_threshold, threshold)
%% This function is to roughly extract the point cloud of the substructure
% Input:
%       Tree                    : Data structure stored origin cells
%       active_pts              : A list of active points
%       road_surface
% %       supstructure            : Data structure of Superstructure
%       struct_threshold
%       threshold
% Output:

% Demo:
%% 
% Tree = OQTR;
% supstructure = SuperStructure;
% active_pts = Active_Pts;
% threshold = THRESHOLD;

%% Extract the leaf nodes
leaf_cell_ids = Node_Leaf(Tree);
% Create a subTree
tempTree = struct('pts', [], 'cell_pts', []);

% Preallocation cell_feature
cells_elev_features = inf(numel(leaf_cell_ids),3);
% Create subTree and calculate cell_height_features
cell_id = 1;
for i = 1:numel(leaf_cell_ids)
    % Retrieve the active points within the cell
    leaf_cell_id = leaf_cell_ids(i);
    leaf_cell_pts_ids = Tree.cell_pts(leaf_cell_id).id;
    
    % Retrieve the remaining points
    leaf_cell_active_pts_ids = leaf_cell_pts_ids(active_pts(leaf_cell_pts_ids));
    
    % Filter the cell
    if numel(leaf_cell_active_pts_ids) >= threshold.min_num_pts
        % Assign the cell to subTree
        tempTree.cell_pts(cell_id).id = leaf_cell_active_pts_ids;
%         subTree.cell_bounds(cell_id,:) = Tree.cell_bounds(leaf_cell_id, :);
        
        % Compute the cell features based on the z-coordinate of the points
        leaf_cell_active_pts_z = Tree.pts(leaf_cell_active_pts_ids,3);
        checked_cell = cell_filtering_pts_elevation(leaf_cell_active_pts_z, struct_threshold);
        cells_elev_features(i,:) = [leaf_cell_id, cell_id, checked_cell];
%             
%         leaf_cell_active_pts_z = sort(Tree.pts(leaf_cell_active_pts_ids,3));
%         cells_elev_features(i,:) = [leaf_cell_id, cell_id,...
%                                     max(diff(leaf_cell_active_pts_z)),...
%                                     abs(leaf_cell_active_pts_z(1) - leaf_cell_active_pts_z(end))];
        % Update cell_id
        cell_id = cell_id + 1;
    end

end
% Remove deactive cells
mask = any(isinf(cells_elev_features),2);
cells_elev_features = cells_elev_features(~mask,:);

% %% Clustering and Filtering the cluster
% % Extract road surface and calculate its feature
% mask = cellfun(@(s) contains(s, 'Road Surface', 'IgnoreCase', true), supstruct(1).Description, 'UniformOutput', false);
% road_surf_id = find(cell2mat(mask));
% road_surf_pts_ids = vertcat(supstruct.Component(road_surf_id).cell.ptc_ids);
% road_surf_pts_xyz = Tree.pts(road_surf_pts_ids,1:3);
% road_mbb = min2DBoundingBox(road_surf_pts_xyz(:,1:2)');
% % road_surf_cent = [mean(road_mbb.vertices,1), 0];
% % road_surf_normal = cross([road_mbb.long_edge_vector, 0], [road_mbb.short_edge_vector,0]);
% road_surf_tangent = [road_mbb.long_edge_vector,0];

%% Clustering cells representing a substructure
% Clustering the cells
substruct_cell_ids = cells_elev_features(logical(cells_elev_features(:,3)),1:2);
substruct_region_cell_info = cell_connectivity_segmentation(Tree, substruct_cell_ids(:,1));

% Check cells in the region      
region_ids = unique(substruct_region_cell_info(:,2));

% Preallocation a score matrix
region_features = zeros(numel(region_ids),5);
region_features(:,1) = region_ids;

% Check progress
for i = 1:numel(region_ids)
    % Retrieve the cells/points within the region by using subTree
    mask = substruct_region_cell_info(:,2) == region_ids(i);
    region_cell_ids = substruct_cell_ids(mask,2);
    region_cell_pts_ids = vertcat(tempTree.cell_pts(region_cell_ids).id);
    region_cell_pts_xyz = Tree.pts(region_cell_pts_ids,1:3);
    region_mbb = min2DBoundingBox(region_cell_pts_xyz(:,1:2)');
    
    % Assign a center of the region
    region_features(i,2:4) = mean(region_cell_pts_xyz,1);
    
    % Establish edge_length and edge_vector of the region_mbb
    edge_length = [region_mbb.long_edge;region_mbb.short_edge];
    edge_vects = [region_mbb.long_edge_vector;region_mbb.short_edge_vector];
    edge_vects(:,3) = road_surface.mbb.long_edge_vector(3); % Project to the same plane of road_surface.mbb.long_edge_vector
    % Compute angles to a traffic direction to find the edge along the
    % longtitudinal direction of the substructure
    cosine_edge_road = cosine_vectors(road_surface.mbb.long_edge_vector, edge_vects);
    [~, min_id] = min(cosine_edge_road);
    longtitudinal_edge_length = edge_length(min_id);
    
    % Filtering component
    if longtitudinal_edge_length >= struct_threshold.traffic_lane_width
        region_features(i,5) = 1;
    end
    
end
region_features = region_features(logical(region_features(:,5)),:);

%% Classify region as the abutment and pier, where the abutments (end points) orther as piers
[~, ~, ~, start_ptc_id, end_ptc_id] = findEndPts(region_features(:,2:4));
abutment_ids = [region_features(start_ptc_id,1), region_features(end_ptc_id,1)];

%% Create data structure 
subTree = struct('Description', {}, 'cell', []);
for i = 1:size(region_features, 1)
    % Retrieve the cell_ids and pts_ids
    mask = substruct_region_cell_info(:,2) == region_features(i,1);
    region_cell_ids = substruct_cell_ids(mask,2);
    
    % Add cells to the substruct
    for j = 1:numel(region_cell_ids)
        subTree(i).cell(j).ids = region_cell_ids(j);
        subTree(i).cell(j).ptc_ids = tempTree.cell_pts(region_cell_ids(j)).id;
    end
    
    % Assign the name
    if any(ismember(abutment_ids, region_features(i, 1)))
        subTree(i).Description = strcat('Abutment', " ", num2str(i));
    else
        subTree(i).Description = strcat('Pier', " ", num2str(i));
    end
end

end % End of the function   