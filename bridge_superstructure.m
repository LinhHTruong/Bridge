function [bridge_tree, cells_region, un_connect_segs_ids, ref_seg_link_list] = bridge_superstructure(Tree, cells_region, bridge_tree, bridge_tree_level, level_name, struct_threshold, threshold, debug)
%% This function is to extract the components surface of the bridge structure
% Input:
%     Tree                      : Data structure stored 2D cell in xy associated with ptc_ids
%     bridge_tree               : Data structure store results
%     bridge_tree_level
%     region                    : Data structure stored points in each segment (plane)
%     comp_name                 : Name of supperstructure
%     struct_threshold          : Structur_threshold relating a minimum size of the structural components
%     threshold                 : Threashold relating to filtering and rocessing point cloud
%     debug                     : Check running time

% Output
%     bridge_tree               : Data structure store results
%     region                    : Data structure stored points in each segment (plane)
%     comp_link_list            : [Nx4] [ref_seg, connect_seg, relationship, side]; 
%                                 relationship: coplanar, parallel, intersection; side = 1 and -1

%% Demo
% Tree = OQTR;
% bridge_tree = Bridge;

% cells_region = Cells_Region;
% level_name = 'top_superstruct';
% struct_threshold = struct_threshold;
% threshold = THRESHOLD;
% debug = true;

%%
Tree = OQTR;
cells_plane = Cells_Plane;
bridge_tree = Bridge;
cells_region = Cells_Region;
level_name = 'bottom_superstruct';
struct_threshold = struct_threshold;
threshold = THRESHOLD;
debug = true;

%% Compute the feature of the region
region_features = cal_region_features(Tree, cells_region);


%% Find the reference component: top: road surface; bottom: bottom surface of girder
if strcmp(level_name, 'top_superstruct')
    % Road surface
    [~, max_id] = max(region_features.bounding_box(:,1));
    ref_seg_ids = region_features.ids(max_id, 1);
%     ref_seg_cell_ids = vertcat(cells_region(ref_seg_ids).cell.id);
    
    % Update bridge tree
    [bridge_tree, ~, cells_region] = bridge_tree_components(bridge_tree, bridge_tree_level, 'Road Surface', 0, ref_seg_ids, cells_region);

elseif strcmp(level_name, 'bottom_superstruct')    
    % Compute a distance from the bottom surfaces to the bridge desk
    bridge_comp_name = (arrayfun(@(x) vertcat(bridge_tree(1).Description(x).ids),1:length(bridge_tree(1).Description),'UniformOutput',false));
    bridge_road_surf_id = find(contains(bridge_comp_name, 'Road Surface'));
    
    % Retrieve the cell within the bridge road surface
    bridge_road_surf_cell_ids = vertcat(bridge_tree(1).Component(bridge_road_surf_id).cell.ids);
    [~, sort_ids] = sort(bridge_road_surf_cell_ids(:,1));
    bridge_road_surf_cell_ids = bridge_road_surf_cell_ids(sort_ids,:);
    mask = ismember(cells_plane.cell_ids(:,[1,2]),bridge_road_surf_cell_ids, 'rows');
    bridge_road_surf_cell_surf = vertcat(cells_plane.peak_info(mask).peaks_features);

    % Compute distance from the bottom surfaces of the superstructure to
    % the road surface by using cell-to-cell method
    no_region = length(cells_region);
    dist_region_road_surface = [];
    for region_id = 1:no_region
        % Retrieve the cells within the region
        cell_ids = vertcat(cells_region(region_id).cell.id); 
        cell_surf_features = vertcat(cells_region(region_id).cell.surface_features); 
        
        % Compute a distance from region to road surface
        dist_cell_road_surface = zeros(size(cell_surf_features,1), 1);
        
        for i = 1:size(cell_ids,1)
            % Search overlap cells
            mask = ismember(bridge_road_surf_cell_ids(:,1), cell_ids(i,1));
            if any(mask)
                % Retrieve a local surface on the bridge road surface: the
                % same cell id
                bridge_road_local_surf = bridge_road_surf_cell_surf(mask,1:6);
            else
                % Retrieve a local surface on the bridge road surface: the
                % neastest cells
                neighbour_cell_id = knnsearch(bridge_road_surf_cell_surf(:,[1,2]),cell_surf_features(i,[1,2]),'K',1);
                bridge_road_local_surf = bridge_road_surf_cell_surf(neighbour_cell_id,1:6);
                clear neighbour_cell_id
            end
            
            % Compute distance: delta_z
%             intersection_pt = cal_line_plane_intersection([cell_surf_features(i,1:3),[0,0,1]], bridge_road_local_surf);
%             proj_ptc = proj_3Dptc_3Dplane(cell_surf_features(i,1:3),bridge_road_local_surf)
            dist_cell_road_surface(i) = abs(cell_surf_features(i,3) - bridge_road_local_surf(3));%abs(dist_3Dpoints_3Dplane(cell_surf_features(i,1:3),bridge_road_local_surf));
            clear mask bridge_road_local_surf
        end
        % Assembly data
%         region_cell_dist = [cell_surf_features(:,[1,2,3]),dist_cell_road_surface, cell_ids];
        dist_cell_road_surface(:,2) = region_id;
        dist_region_road_surface = [dist_region_road_surface;region_cell_dist];
        clear cell_ids cell_surf_features dist_cell_road_surface
    end 


    
    
    
    
    
    
    
    
    
else 
    
end

%% Find the region connect to the reference surface
if strcmp(level_name, 'top_superstruct')
    mask = (region_features.center(:,3) >= region_features.center(ref_seg_ids,3))&...
           (region_features.bounding_box(:,2) >= 0.5*struct_threshold.bridge_min_length);
    connect_segs_ids = region_features.ids(mask,1);
    connect_segs_ids = setdiff(connect_segs_ids,ref_seg_ids); 
    un_connect_segs_ids = region_features.ids(~mask,1);
    comp_link_list = connected_2D_components(Tree, cells_region, region_features,...
                            ref_seg_ids, connect_segs_ids, struct_threshold, threshold);
                
%     % component_link_list = [master id, slave id, type connection[-1, 0, 1], side connection[-1,1]               
%     % Update un_components_regions_ids                    
%     if isempty(un_connect_segs_ids)
%         un_connect_segs_ids = setdiff(connect_segs_ids, comp_link_list(:,[1,2]));
%     else
%         un_connect_segs_ids = union(un_connect_segs_ids, setdiff(connect_segs_ids, comp_link_list(:,[1,2])));
%     end
%     un_connect_segs_ids = setdiff(un_connect_segs_ids, ref_seg_cell_ids); 

    % Check if child on the same side and connect together
    flag = true(size(comp_link_list,1),1);
    for i = 1:size(comp_link_list,1)
        connect_comp_ids = comp_link_list(i,:);
        mask = ismember(comp_link_list(:,[1,2]),connect_comp_ids(1:2), 'rows');
        remain_comp_link_list = comp_link_list(~mask,:);
        mask = ismember(remain_comp_link_list(:,2), connect_comp_ids(:));
        parent_id = (remain_comp_link_list(mask, 1));
        side_direction = (remain_comp_link_list(mask, 3));
        if (numel(parent_id) == 2)&&(numel(unique(parent_id)) == 1)&&(numel(side_direction) == 2)&&(numel(unique(side_direction) == 1))
            flag(i) = false;
        end

    end
    comp_link_list = comp_link_list(flag,:);
    clear i mask connect_component_ids remain_component_link_list parent_id side_direction flag

    % Assign comp_link_list to the bridge component
    mask = ismember(comp_link_list(:,1), ref_seg_ids);
    ref_seg_link_list = comp_link_list(mask,:);
    non_ref_seg_link_list = comp_link_list(~mask,:);
    
    % Update the ref_seg_linked_list
    unique_link = unique(ref_seg_link_list(:,4));
    for i = 1:numel(unique_link)
        % Retrieve segments connected to the ref_seg
        mask = ref_seg_link_list(:,4) == unique_link(i);
        conected_seg_ids = ref_seg_link_list(mask,2);
        [bridge_tree, ~, cells_region] = bridge_tree_components(bridge_tree, bridge_tree_level, 'Foot Path', ref_seg_ids, conected_seg_ids, cells_region);
    end
    
    % Update non-ref_seg_link_list
    unique_link = unique(non_ref_seg_link_list(:,4));
    for i = 1:numel(unique_link)
        % Retrieve segments connected to the ref_seg
        mask = non_ref_seg_link_list(:,4) == unique_link(i);
        conected_seg_ids = non_ref_seg_link_list(mask,2);
        [bridge_tree, ~, cells_region] = bridge_tree_components(bridge_tree, bridge_tree_level, 'Mis', ref_seg_ids, conected_seg_ids, cells_region);
    end
    

elseif strcmp(level_name, 'immediate_superstruct')
    
    
elseif strcmp(level_name, 'bottom_superstruct')
    
else
    % Error
end





