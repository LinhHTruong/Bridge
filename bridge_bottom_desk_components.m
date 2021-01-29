function [bridge, un_components_regions_ids] = bridge_bottom_desk_components(Tree, Cells, Cell_Region, Region_Features, ...
                                                bridge, bridge_tree_level, threshold, bridge_threshold)
%% This function is to extract the points cloud of the components connected to bridge surface
% - Bridge surface: the largest region
% - Footpath: 1 or 2 or 0 beside the bridge surface: (i) location: above
% the bridge surface; (ii) width >= 1.0m; (iii) length > 0.25 bridge
% surface

% - Parapet: Beside the footpath
% - guard rail between the footpath and bridge desk

% Demo:
% Tree = OQTR;
% Cells = Plane_Cells;
% Cell_Region = Region;
% bridge = Bridge;
% bridge_tree_level = 2;

%% Extract the bridge desk
if isempty(bridge(1).Component)
    fprintf('Bridge component structure is empty');
    return
else
    % Extract points and compute features of the bridge desk
    num_cell = size(vertcat(bridge(1).Component(1).cell.ids),1);
    bridge_deck_cells_center = zeros(num_cell,3);
    for j = 1:num_cell
        cell_peak_ptc_ids = bridge(1).Component(1).cell(j).ptc_ids;
        cell_peak_ptc_xyz = Tree.pts(cell_peak_ptc_ids,1:3);
        bridge_deck_cells_center(j,:) = mean(cell_peak_ptc_xyz,1);
        clear j cell_peak_ptc_ids cell_peak_ptc_xyz
    end
        
    bridge_deck_center = mean(bridge_deck_cells_center,1);
    [bridge_deck_eigen_vectors, ~, ~] = eigenspace(bridge_deck_cells_center, 0);
    bridge_deck_normal = bridge_deck_eigen_vectors(1,:);
%     
%     bridge_deck_cell_peak_ptc_ids = vertcat(bridge(1).Component(1).cell.ptc_ids);
%     bridge_deck_cell_peak_ptc_xyz = TREE.pts(bridge_deck_cell_peak_ptc_ids,1:3);
%     bridge_deck_center = mean(bridge_deck_cell_peak_ptc_xyz,1);
%     [bridge_deck_eigen_vectors, ~, ~] = eigenspace(bridge_deck_cell_peak_ptc_xyz, 0);
%     bridge_deck_normal = bridge_deck_eigen_vectors(1,:);
%     
end
%% Compute orthogonal distance from the cells in the region to the bridge desk
no_region = length(Cell_Region);
regions_cells_peaks_distance = [];
for region_id = 1:no_region
    % Extract the points in the peaks and compute distances from the peak
    % to the bridge desk
    cell_peak_ids = vertcat(Cell_Region(region_id).cell.id); 
    cell_peak_surface_features = vertcat(Cell_Region(region_id).cell.surface_features); 
    cell_peak_ids(:,3) = region_id;
    cell_peak_ids(:,4) = abs(dist_3Dpoints_3Dplane(cell_peak_surface_features(:,1:3), [bridge_deck_center, bridge_deck_normal]));
    % Cummulative region
    regions_cells_peaks_distance = [regions_cells_peaks_distance;cell_peak_ids];
end 
clear region_id cell_peak_ids cell_peak_surface_features
%% Determine the regions are the bottom surface of the bridge desk by using kde: the distance is mostly constant: DBscan is another option for curvature
[fi,zi,~] = ksdensity(regions_cells_peaks_distance(:,end),'npoints',100,'bandwidth',threshold.distance,'Kernel','epanechnikov');
peak_shape = peak_shape_width(zi, fi);

% Determine cell_peak having the distance within the largest peak
[~, id] = max(fi);
largest_peak_zi = zi(id);
mask = peak_shape(:,1) == largest_peak_zi;
largest_peak_width = peak_shape(mask,2);

mask = (largest_peak_zi - largest_peak_width/2 <= regions_cells_peaks_distance(:,end))&...
       (regions_cells_peaks_distance(:,end) <= largest_peak_zi + largest_peak_width/2);
   
largest_peak_regions_cell_peaks_ids = regions_cells_peaks_distance(mask,1:3);
clear fi zi peak_shape id mask largest_peak_zi largest_peak_width
%% Decision which regions are the bottom of the bridge desk
region_ids = unique(largest_peak_regions_cell_peaks_ids(:,end));
flag_region = false(numel(region_ids), 1);
for i = 1:numel(region_ids)
    % The region is the bottom surface if 80% of the cell peaks are within
    % the largest peak in term of the distance
    % From the distance
    mask = largest_peak_regions_cell_peaks_ids(:,end) == region_ids(i);
    count_detected_cell_peaks = sum(mask);
    
    % From the region
    cell_peak_ids = vertcat(Cell_Region(region_ids(i)).cell.id); 
    if count_detected_cell_peaks >= 0.75*size(cell_peak_ids,1)
        flag_region(i) = true;
    end
end
bridge_bottom_deck_regions_ids = region_ids(flag_region);   
clear region_ids flag_region i mask count_detected_cell_peaks cell_peak_ids
%% Establish the root node
bridge_bottom_deck_surface_feature = zeros(numel(bridge_bottom_deck_regions_ids),6);
for i=1:numel(bridge_bottom_deck_regions_ids)
    
    [bridge,~, Cell_Region] = bridge_tree_components(bridge, bridge_tree_level, 0, bridge_bottom_deck_regions_ids(i), Cell_Region);
    cell_peak_ptc_ids = vertcat(Cell_Region(bridge_bottom_deck_regions_ids(i)).cell.ptc_ids); 
    cell_peak_ptc_xyz = Tree.pts(cell_peak_ptc_ids,1:3);
    bridge_bottom_deck_surface_feature(i,:) = [mean(cell_peak_ptc_xyz, 1), eigenspace(cell_peak_ptc_xyz,1)];
end
clear i cell_peak_ptc_ids cell_peak_ptc_xyz
%% Find connected regions of the bridge desk: Footpath
% mask = (Region_Features.ids ~= bridge_surface_region_ids);
mask = Region_Features.center(:,3) > max(Region_Features.center(bridge_bottom_deck_regions_ids,3));
candidate_components_segments_ids = Region_Features.ids(mask,1);
un_components_regions_ids = Region_Features.ids(~mask,1);

component_link_list = connected_2D_components(Tree, Cell_Region, Region_Features,...
                        bridge_bottom_deck_regions_ids, candidate_components_segments_ids, threshold, bridge_threshold);

% Update un_components_regions_ids                    
if isempty(un_components_regions_ids)
    un_components_regions_ids = setdiff(candidate_components_segments_ids, component_link_list(:,[1,2]));
else
    un_components_regions_ids = union(un_components_regions_ids, setdiff(candidate_components_segments_ids, component_link_list(:,[1,2])));
end
un_components_regions_ids = setdiff(un_components_regions_ids,bridge_bottom_deck_regions_ids); 
%% Update bridge tree Update bridge data structure:
for i = 1:size(component_link_list,1)
    [bridge, ~, Cell_Region] = bridge_tree_components(bridge, bridge_tree_level, component_link_list(i,1), component_link_list(i,2), Cell_Region);
end
    
%% Extract the points between the bottom and a top deck surfaces
% Extract cells on leaf node components: leaf node component is on the boundary of the bottom surface
mask = (bridge(bridge_tree_level).Link_List(:,1) ~= 0);
leaf_nodes_component_ids = setdiff(bridge(bridge_tree_level).Link_List(mask,2), bridge(bridge_tree_level).Link_List(mask,1));
leaf_nodes_component_cells_peaks_info = (arrayfun(@(x) vertcat(Cell_Region(x).cell.id),leaf_nodes_component_ids,'UniformOutput',false));
leaf_nodes_component_cells_peaks_info = cell2mat(leaf_nodes_component_cells_peaks_info);
leaf_nodes_component_cells_ids = unique(leaf_nodes_component_cells_peaks_info(:,1));
clear leaf_nodes_component_cells_peaks_info

% Find the cell connecting to the leaf_node components: cells do not
% assigned to the components and connected to the components of the bottom surface of the bridge deck

% Extract cell_peaks of the bridge components in the bottom suface
num_component = length(bridge(bridge_tree_level).Component);
bridge_bottom_component_cells_peaks_ids = [];
for i = 1:num_component
    % Retrieve cell peak ids in the component
    cell_peak_ids = vertcat(bridge(bridge_tree_level).Component(i).cell.ids);
    cell_peak_ids(:,3) = i;
    bridge_bottom_component_cells_peaks_ids = [bridge_bottom_component_cells_peaks_ids;cell_peak_ids];
end
clear cell_peak_ids

% Extract cell_peaks have not yet assigned to the components 
non_assign_cell_peaks_ids = Cells.cell_ids(Cells.cell_ids(:, end) == 1, 1:2); %Non-assign before this step
non_component_cell_ids = setdiff(non_assign_cell_peaks_ids(:,[1,2]), bridge_bottom_component_cells_peaks_ids(:,[1,2]),'rows', 'stable');

% Searching cell connected to the boundary cells
non_component_connect_cell_ids = Query_Neighbour_Cells(Tree, leaf_nodes_component_cells_ids,non_component_cell_ids);

% 
% 
% flag = false(numel(non_component_cell_ids),1);
% for i = 1:numel(leaf_nodes_component_cells_ids)
%     neighbour_cell_ids = Window_Query_Neighbour_Cell(Tree, leaf_nodes_component_cells_ids(i), non_component_cell_ids, bridge_threshold.search_scale);
%     mask = ismember(non_component_cell_ids,neighbour_cell_ids);
%     flag(mask) = true;
% end
% non_component_connect_cell_ids = non_component_cell_ids(flag);
clear non_assign_cell_peaks_ids non_component_cell_ids neighbour_cell_ids mask flag
%% Find the points between the bottom and top surfaces
% Extract the cell_peaks assigned to a top surface of the bridge
num_component = length(bridge(1).Component);
bridge_top_deck_cell_peak_ids = [];
for i = 1:num_component
    % Retrieve cell peak ids in the component
    cell_peak_ids = vertcat(bridge(1).Component(i).cell.ids);
    if size(cell_peak_ids, 2) == 2
        cell_peak_ids(:,3) = i;
        bridge_top_deck_cell_peak_ids = [bridge_top_deck_cell_peak_ids;cell_peak_ids];
    end
end
clear i num_component cell_peak_ids

% Extract the cell_peaks assigned to a bottom surface of the bridge
bridge_bottom_deck_cell_ids = union(unique(bridge_bottom_component_cells_peaks_ids(:,1)),non_component_connect_cell_ids);

%% Find the points between two surface
connect_component_cells = struct('cell',[]);
count = 0;
for i=1:numel(bridge_bottom_deck_cell_ids)
    
    % Searching the cell_peak of the top deck having a cell id corresponding to the the bottom deck
    mask = ismember(bridge_top_deck_cell_peak_ids(:,1), bridge_bottom_deck_cell_ids(i));
    top_deck_cell_peak_component_ids = bridge_top_deck_cell_peak_ids(mask,:);
    
    if isempty(top_deck_cell_peak_component_ids)
        % This cell is not available in the top surface of the bridge deck
        % due to missing data or over-hange -> searching neighbour
        neighbour_cell_ids = Query_Neighbour_Cells(Tree, bridge_bottom_deck_cell_ids(i), bridge_top_deck_cell_peak_ids(:,1));
        mask = ismember(bridge_top_deck_cell_peak_ids(:,1), neighbour_cell_ids);
        top_deck_cell_peak_component_ids = bridge_top_deck_cell_peak_ids(mask,:);
        clear neighbour_cell_ids mask
    end
        
    % Retrieve the points within top_cell_peak_component_ids and estimate a
    % local surface of a top deck
    component_ids = unique(top_deck_cell_peak_component_ids(:,3));
    top_deck_cell_peak_ptc_ids = [];
    for j = 1:numel(component_ids)
        component_id = component_ids(j);
        cell_peak_component_id = top_deck_cell_peak_component_ids(top_deck_cell_peak_component_ids(:,3) == component_id,1:2);
        mask = ismember(vertcat(bridge(1).Component(component_id).cell.ids), cell_peak_component_id,'rows');
        top_deck_cell_peak_ptc_ids = [top_deck_cell_peak_ptc_ids;vertcat(bridge(1).Component(component_ids(j)).cell(mask).ptc_ids)];
    end
    clear j component_ids component_id cell_peak_component_id mask 
    
    % A local surface of the bridge top deck
    if numel(top_deck_cell_peak_ptc_ids) < threshold.min_num_pts%isempty(top_deck_cell_ptc_ids)
        bridge_top_deck_local_surface = [bridge_deck_center, bridge_deck_normal];
    else
        cell_peak_component_ptc_xyz = Tree.pts(top_deck_cell_peak_ptc_ids,1:3);
        bridge_top_deck_local_surface = [mean(cell_peak_component_ptc_xyz,1), eigenspace(cell_peak_component_ptc_xyz,1)];
        % Add a condition if the the cell is sigularity which contains the
        % points on the vertical surface
        if abs(cosine_vectors(bridge_top_deck_local_surface(:,4:6), [0., 0., 1.])) <= cos(deg2rad(threshold.max_hor_plane_angle))
            bridge_top_deck_local_surface = [bridge_deck_center, bridge_deck_normal];
        end  
        clear cell_peak_component_ptc_xyz
    end
            
    % Find the local surface of the bridge bottom deck
    mask = ismember(bridge_bottom_component_cells_peaks_ids(:,1), bridge_bottom_deck_cell_ids(i)); % This condition is to avoid the cell is connected cells
    if any(mask)
        % The cell belong to the components
        bottom_deck_cell_peak_component_ids = bridge_bottom_component_cells_peaks_ids(mask,:);
    else
        % Searching neighbour of the current cell
%         neighbour_cell_ids = Query_Neighbour_Cells(Tree, bridge_bottom_deck_cell_ids(i), bridge_bottom_component_cells_peaks_ids(:,1));
        neighbour_cell_ids = Window_Query_Neighbour_Cell(Tree, bridge_bottom_deck_cell_ids(i), bridge_bottom_component_cells_peaks_ids(:,1), bridge_threshold.search_scale);
        mask = ismember(bridge_top_deck_cell_peak_ids(:,1), neighbour_cell_ids);
        bottom_deck_cell_peak_component_ids = bridge_top_deck_cell_peak_ids(mask,:);
        clear neighbour_cell_ids mask
    end
    
    % Estimate the local surface of the bridge bottom deck
    if ~isempty(bottom_deck_cell_peak_component_ids)
        component_ids = unique(bottom_deck_cell_peak_component_ids(:,3));
        bottom_deck_cell_peak_ptc_ids = [];
        for j = 1:numel(component_ids)
            component_id = component_ids(j);
            cell_peak_component_id = bottom_deck_cell_peak_component_ids(bottom_deck_cell_peak_component_ids(:,3) == component_id,1:2);
            mask = ismember(vertcat(bridge(bridge_tree_level).Component(component_id).cell.ids), cell_peak_component_id,'rows');
            bottom_deck_cell_peak_ptc_ids = [bottom_deck_cell_peak_ptc_ids;vertcat(bridge(bridge_tree_level).Component(component_ids(j)).cell(mask).ptc_ids)];
        end
        clear component_ids component_id bottom_deck_cell_peak_component_ids j 
        
        % A local surface
        if numel(bottom_deck_cell_peak_ptc_ids) < 20 %isempty(bottom_deck_cell_ptc_ids)
            [~, mask] = min(bridge_bottom_deck_surface_feature(:,3));
            bridge_bottom_deck_local_surface = bridge_bottom_deck_surface_feature(mask,:);
        else
            cell_peak_component_ptc_xyz = Tree.pts(bottom_deck_cell_peak_ptc_ids,1:3);
            bridge_bottom_deck_local_surface = [mean(cell_peak_component_ptc_xyz,1), eigenspace(cell_peak_component_ptc_xyz,1)];
            clear cell_peak_component_ptc_xyz
        end
        
    else
        [~, mask] = min(bridge_bottom_deck_surface_feature(:,3));
        bridge_bottom_deck_local_surface = bridge_bottom_deck_surface_feature(mask,:);
        bottom_deck_cell_peak_ptc_ids = [];
    end

    % Retrieve the points in the original cells and remove points in the bottom and top deck
    original_cell_ptc_ids = Tree.cell_pts(bridge_bottom_deck_cell_ids(i)).id;
    remain_cell_ptc_ids = setdiff(original_cell_ptc_ids, union(top_deck_cell_peak_ptc_ids, bottom_deck_cell_peak_ptc_ids));
    remain_cell_ptc_xyz = Tree.pts(remain_cell_ptc_ids,1:3);
    clear bottom_deck_cell_peak_ptc_ids top_deck_cell_peak_ptc_ids
    
    % Compute distance from these points to the local surface of the bridge top deck
    dist_ptc_top_surface = dist_3Dpoints_3Dplane(remain_cell_ptc_xyz, bridge_top_deck_local_surface);
    dist_bottom_top_surface = dist_3Dpoints_3Dplane([bridge_top_deck_local_surface(1:2),bridge_bottom_deck_local_surface(3)], bridge_top_deck_local_surface);
            
    % Filter the points between two surfaces + tolerance
    threshold_dist_range = [+sign(dist_bottom_top_surface)*threshold.residual, -sign(dist_bottom_top_surface)*threshold.residual + dist_bottom_top_surface];
    mask = (min(threshold_dist_range) <= dist_ptc_top_surface)&(dist_ptc_top_surface <= max(threshold_dist_range));
    bottom_top_deck_ptc_ids = remain_cell_ptc_ids(mask);
    
    if numel(bottom_top_deck_ptc_ids) >= 20
        bottom_top_deck_ptc_xyz = remain_cell_ptc_xyz(mask,1:3);
%         connect_bottom_top_deck_height = max(connect_bottom_top_deck_ptc_xyz(:,3)) - min(connect_bottom_top_deck_ptc_xyz(:,3));

        % Assign to the connect_component_region
        connect_component_cells.cell(count + 1).id = bridge_bottom_deck_cell_ids(i);
        connect_component_cells.cell(count + 1).center = mean(bottom_top_deck_ptc_xyz,1);
        connect_component_cells.cell(count + 1).ptc_ids = bottom_top_deck_ptc_ids;
%         connect_component_cells.cell(count + 1).height = connect_bottom_top_deck_height;
        count = count + 1;
    end 
    clear bridge_top_deck_local_surface bridge_bottom_deck_local_surface threshold_dist_range dist_ptc_top_surface dist_bottom_top_surface
    clear bottom_top_deck_ptc_ids bottom_top_deck_ptc_xyz
end

%% Apply clustering technique to remove incorrect cells
cell_ids = vertcat(connect_component_cells.cell.id);
cell_center = vertcat(connect_component_cells.cell.center);
% cell_normal = vertcat(connect_component_cells.cell.normal);
% cell_height = vertcat(connect_component_cells.cell.height);

% Compute features of the cells for outlier removal, segmentation: distance from cell_proj-center to the central line of the bridge deck
proj_cell_center = proj_3Dpoints_3Dplane(cell_center, [bridge_deck_center,bridge_deck_normal]);
dist_cell_central_line = dist_3Dpoints_3Dline(proj_cell_center, [bridge_deck_center, bridge_deck_eigen_vectors(3,:)], bridge_deck_eigen_vectors(2,:));
dist_cell_bridge_deck = dist_3Dpoints_3Dplane(cell_center, [bridge_deck_center,bridge_deck_normal]);

% Use DB scan for distance from the cell to the bridge surface
dist_cell_bridge_deck_epsilon = 2.0*threshold.sampling_step; % Parallel -> cells in the same component having the same distance
dist_cell_bridge_deck_minpts = 0.5*bridge_threshold.min_no_cells; %
db_dist_cell_bridge_deck = clustering_dbscan(dist_cell_bridge_deck, dist_cell_bridge_deck_epsilon, dist_cell_bridge_deck_minpts, -1);

% Use DB scan for sub-group
segment_cell = zeros(numel(cell_ids),2);
segment_cell(:,1) = cell_ids;
num_cluster = max(db_dist_cell_bridge_deck);

count_segment = 1;
for j = 1:num_cluster
    mask = db_dist_cell_bridge_deck == j;
    cluster_cell_ids = cell_ids(mask);
    cluster_cell_center = cell_center(mask,:);
    cluster_dist_cell_central_line = dist_cell_central_line(mask);
    
    % Dbscan for dist_cell_2_central_line
    db_dist_cell_central_line = clustering_dbscan(cluster_dist_cell_central_line, dist_cell_bridge_deck_epsilon, dist_cell_bridge_deck_minpts, -1);

    % Check sub-cluster
    num_subcluster = max(db_dist_cell_central_line);
    
    for k = 1:num_subcluster
        % Extract the cell in the sub-cluster
        mask = db_dist_cell_central_line == k;
        subcluster_cell_ids = cluster_cell_ids(mask,:);
        subcluster_cell_center = cluster_cell_center(mask,:);
        
        % Compute the length of the sub-cluster
        [~, ~, endPoints, ~, ~] = findEndPts(subcluster_cell_center);
        subcluster_length = norm(endPoints(1:3) - endPoints(4:6));
        
        if subcluster_length >= dist_cell_bridge_deck_minpts*threshold.cell_size
            mask = ismember(segment_cell(:,1), subcluster_cell_ids);
            segment_cell(mask,2) = count_segment;
            count_segment = count_segment + 1;
        end
    end
end

%% Assign the cluster for 
component_region = struct('cell',[], 'status',[]);
mask = segment_cell(:,2) > 0;
segment_cell = segment_cell(mask,:);
for j = 1:max(segment_cell(:,2))
    mask = segment_cell(:,2) == j;
    segment_cell_ids = segment_cell(mask,1);
    
    % Update struct
    component_region(j).status = 1;
    for count = 1:numel(segment_cell_ids)
        mask = ismember(cell_ids, segment_cell_ids(count));
        component_region(j).cell(count).id = connect_component_cells.cell(mask).id;
        component_region(j).cell(count).ptc_ids = connect_component_cells.cell(mask).ptc_ids;
    end
end

%% Update bridge data structure: child of child
bridge(bridge_tree_level).Link_List(:,end) = 0;
for k = 1:length(component_region)
    [bridge, ~, component_region] = bridge_tree_components(bridge, bridge_tree_level, bridge_bottom_deck_regions_ids(1), k, component_region);
end
bridge(bridge_tree_level).Link_List(:,end) = [];
