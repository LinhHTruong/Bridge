
function [bridge, un_components_segments_ids] = bridge_desk_components(Tree, Cell_Region, Region_Features, bridge, bridge_tree_level, threshold, bridge_threshold)
%% This function is to extract the points cloud of the components connected to bridge surface
% - Bridge surface: the largest region
% - Footpath: 1 or 2 or 0 beside the bridge surface: (i) location: above
% the bridge surface; (ii) width >= 1.0m; (iii) length > 0.25 bridge
% surface

% - Parapet: Beside the footpath
% - guard rail between the footpath and bridge desk

% Demo:
% Tree = OQTR;
% bridge = Bridge
% bridge_tree_level = 1;
% Cell_Region = Region;
% threshold = THRESHOLD

%% Extract leaf node
leaf_cell_ids = Node_Leaf(Tree);
leaf_cell_ids = leaf_cell_ids(Tree.cell_props(leaf_cell_ids) == 1);

%% Find the largest plane in term of area: the bridge desk
[~, max_id] = max(Region_Features.bounding_box(:,1));
bridge_deck_segment_ids = Region_Features.ids(max_id, 1);
bridge_deck_cell_ids = vertcat(Cell_Region(bridge_deck_segment_ids).cell.id);

% Update bridge data structure
[bridge, ~, Cell_Region] = bridge_tree_components(bridge, bridge_tree_level, 0, bridge_deck_segment_ids, Cell_Region);

%% Find connected regions of the bridge desk: Footpath
% mask = (Region_Features.ids ~= bridge_surface_region_ids);
mask = Region_Features.center(:,3) > Region_Features.center(bridge_deck_segment_ids,3);
candidate_components_segments_ids = Region_Features.ids(mask,1);
un_components_segments_ids = Region_Features.ids(~mask,1);
component_link_list = connected_2D_components(Tree, Cell_Region, Region_Features,...
                        bridge_deck_segment_ids, candidate_components_segments_ids, threshold, bridge_threshold);
% component_link_list = [master id, slave id, type connection[-1, 0, 1], side connection[-1,1]               
% Update un_components_regions_ids                    
if isempty(un_components_segments_ids)
    un_components_segments_ids = setdiff(candidate_components_segments_ids, component_link_list(:,[1,2]));
else
    un_components_segments_ids = union(un_components_segments_ids, setdiff(candidate_components_segments_ids, component_link_list(:,[1,2])));
end
un_components_segments_ids = setdiff(un_components_segments_ids, bridge_deck_segment_ids); 

% % Update bridge data structure: See how to link 
% 
% bridge(bridge_tree_level).Link_List(:,end) = 0;
% for i = 1:size(component_link_list)
%     [bridge, ~, Cell_Region] = bridge_tree_components(bridge, bridge_tree_level, component_link_list(i,1), component_link_list(i,2), Cell_Region);
%       
% end   
% Check if child on the same side and connect together
flag = true(size(component_link_list,1),1);
for i = 1:size(component_link_list,1)
    connect_component_ids = component_link_list(i,:);
    mask = ismember(component_link_list(:,[1,2]),connect_component_ids(1:2), 'rows');
    remain_component_link_list = component_link_list(~mask,:);
    mask = ismember(remain_component_link_list(:,2), connect_component_ids(:));
    parent_id = (remain_component_link_list(mask, 1));
    side_direction = (remain_component_link_list(mask, 3));
    if (numel(parent_id) == 2) &(numel(unique(parent_id)) == 1)&(numel(side_direction) == 2)&(numel(unique(side_direction) == 1))
        flag(i) = false;
    end
        
end
component_link_list = component_link_list(flag,:);
clear i mask connect_component_ids remain_component_link_list parent_id side_direction flag

%% Extract the point between two parallel components
for i = 1:size(component_link_list)
    if component_link_list(i,3) == -1 % parallel component
       
        % Preallocation
        connect_component_region = struct('cell',[], 'status', 1);
        count = 0;

        % Retrieve the cells of the master/parent component and its bound cells
        parent_id = component_link_list(i,1);
        parent_cell_ids = vertcat(Cell_Region(parent_id).cell.id);
        [parent_bound_cell_ids, ~] = boundary_cells(Tree, parent_cell_ids(:,1));
        parent_non_bound_cell_ids = setdiff(parent_cell_ids(:,1),parent_bound_cell_ids); 
        
        % Retrieve the cells of the slave/child component and its bound cells
        child_id = component_link_list(i,2);
        child_cell_ids = vertcat(Cell_Region(child_id).cell.id);
        [child_bound_cell_ids, ~] = boundary_cells(Tree, child_cell_ids(:,1));
        child_non_bound_cell_ids = setdiff(child_cell_ids(:,1),child_bound_cell_ids);
        
        % Find neighbour cells of the bound cells of the parent: bound
        % cells of the childand others; similar apply to the child
        reference_cell_ids = setdiff(leaf_cell_ids, parent_cell_ids(:,1));
        parent_bound_neighbour_cell_ids = Query_Neighbour_Cells(Tree, parent_bound_cell_ids, reference_cell_ids);
        parent_bound_neighbour_cell_ids = setdiff(parent_bound_neighbour_cell_ids,parent_non_bound_cell_ids); 
        
        reference_cell_ids = setdiff(leaf_cell_ids, child_cell_ids(:,1));
        child_bound_neighbour_cell_ids = Query_Neighbour_Cells(Tree, child_bound_cell_ids,reference_cell_ids);
        child_bound_neighbour_cell_ids = setdiff(child_bound_neighbour_cell_ids,child_non_bound_cell_ids); 
        
        % Extract the points on the boundary of two components
        bound_cell_ids = intersect(parent_bound_neighbour_cell_ids, child_bound_neighbour_cell_ids);
        bound_cell_ids = (union(bound_cell_ids, intersect(parent_bound_cell_ids,child_bound_cell_ids)));
        
        % Filter the points
        if ~isempty(bound_cell_ids)
            for j = 1:numel(bound_cell_ids)

                check_cell_id = bound_cell_ids(j);
                % Estimate a local surface at the parent component
                mask = ismember(parent_cell_ids(:,1),check_cell_id);
                if any(mask)
                    parent_check_cell_ptc_ids = Cell_Region(parent_id).cell(mask).ptc_ids;
                    parent_local_surface = Cell_Region(parent_id).cell(mask).surface_features;
                else
                    % Searching neighbour cells and estimate a local
                    % surface based on the points within these cells
                    parent_check_cell_ptc_ids = [];
                    parent_neighbour_check_cell_ids = Query_Neighbour_Cells(Tree, check_cell_id, parent_cell_ids);
                    mask = ismember(parent_cell_ids(:,1),parent_neighbour_check_cell_ids);
                    parent_neighbour_check_cell_ptc_ids = vertcat(Cell_Region(parent_id).cell(mask).ptc_ids);
                    parent_neighbour_check_cell_ptc_xyz = Tree.pts(parent_neighbour_check_cell_ptc_ids,1:3);
                    parent_local_surface = [mean(parent_neighbour_check_cell_ptc_xyz, 1), eigenspace(parent_neighbour_check_cell_ptc_xyz,1)];
                    clear mask parent_neighbour_check_cell_ids parent_neighbour_check_cell_ptc_ids parent_neighbour_check_cell_ptc_xyz
                end
                
                % Estimate a local surface at the child component
                mask = ismember(child_cell_ids(:,1),check_cell_id);
                if any(mask)
                    child_check_cell_ptc_ids = Cell_Region(child_id).cell(mask).ptc_ids;
                    child_local_surface = Cell_Region(child_id).cell(mask).surface_features;
                else
                    % Searching neighbour cells and estimate a local
                    % surface based on the points within these cells
                    child_check_cell_ptc_ids = [];
                    child_neighbour_check_cell_ids = Query_Neighbour_Cells(Tree, check_cell_id, child_cell_ids);
                    mask = ismember(child_cell_ids(:,1),child_neighbour_check_cell_ids);
                    child_neighbour_check_cell_ptc_ids = vertcat(Cell_Region(child_id).cell(mask).ptc_ids);
                    child_neighbour_check_cell_ptc_xyz = Tree.pts(child_neighbour_check_cell_ptc_ids,1:3);
                    child_local_surface = [mean(child_neighbour_check_cell_ptc_xyz, 1), eigenspace(child_neighbour_check_cell_ptc_xyz,1)];
                    clear mask child_neighbour_check_cell_ids child_neighbour_check_cell_ptc_ids child_neighbour_check_cell_ptc_xyz
                end
                
                % Retrieve the points within the check-cell_id and remove
                % any points within a cell both/either parent and/or child component
                check_cell_all_ptc_ids = Tree.cell_pts(check_cell_id).id;
                if ~isempty(parent_check_cell_ptc_ids)
                    check_cell_all_ptc_ids = setdiff(check_cell_all_ptc_ids, parent_check_cell_ptc_ids);
                end
                if ~isempty(child_check_cell_ptc_ids)
                    check_cell_all_ptc_ids = setdiff(check_cell_all_ptc_ids, child_check_cell_ptc_ids);
                end
                check_cell_all_ptc_xyz = Tree.pts(check_cell_all_ptc_ids, 1:3);
                
                % Extract the points between the two surfaces
                dist_ptc_parent_surface = dist_3Dpoints_3Dplane(check_cell_all_ptc_xyz, parent_local_surface);
%                 dist_parent_child_surface = dist_3Dpoints_3Dplane(parent_local_surface(1:3), child_local_surface);
                dist_child_parent_surface = dist_3Dpoints_3Dplane(child_local_surface(1:3), parent_local_surface);
                
                threshold_dist_range = [+sign(dist_child_parent_surface)*threshold.residual, -sign(dist_child_parent_surface)*threshold.residual + dist_child_parent_surface];
                mask = (min(threshold_dist_range) <= dist_ptc_parent_surface)&(dist_ptc_parent_surface <= max(threshold_dist_range));
                connect_parent_child_cell_ptc_ids = check_cell_all_ptc_ids(mask);
                
                % Update components
                if numel(connect_parent_child_cell_ptc_ids) >= threshold.min_num_pts           
                    % Assign to the connect_component_region
                    connect_component_region.cell(count + 1).id = check_cell_id;
                    connect_component_region.cell(count + 1).center = mean(check_cell_all_ptc_xyz(mask,1:3),1);
                    connect_component_region.cell(count + 1).ptc_ids = connect_parent_child_cell_ptc_ids;
                    count = count + 1;
                end
            end
        end
        
        % Update bridge component
        bridge(bridge_tree_level).Link_List(:,end) = 0;
        if ~isempty(connect_component_region.cell)
%             [bridge, ~] = bridge_tree_components(bridge, bridge_tree_level, parent_id, 1, connect_component_region);
            [bridge, ~, connect_component_region] = bridge_tree_components(bridge, bridge_tree_level, parent_id, 1, connect_component_region);

        end      
        
%     elseif component_link_list(i,3) == 1 % interaction component -> filtering points outside the boundary
        %Will implement in future
    end
end

%% Determine parapets and safeguard rails
reference_cell_ids = setdiff(leaf_cell_ids, bridge_deck_cell_ids);

parent_ids = unique(component_link_list(:,1));

for i = 1:numel(parent_ids)
    %%
    parent_id = parent_ids(i);
    mask = component_link_list(:,1) == parent_id;
    sub_component_link_list = component_link_list(mask,:);
    sub_components_sign = unique(sub_component_link_list(:,end));

    for j = 1:numel(sub_components_sign)
        % Retrieve region on one side
        mask = sub_component_link_list(:,end) == sub_components_sign(j);
        child_component_ids = unique(sub_component_link_list(mask,2));
        % Combinate all child components
        components_struct = merge_structs(Cell_Region, child_component_ids);
        components_struct.status = 1;
        
        % Update bridge tree
        bridge(bridge_tree_level).Link_List(:,end) = 0;
        [bridge, new_parent_id, components_struct] = bridge_tree_components(bridge, bridge_tree_level, parent_id, 1, components_struct);
    %      
        % Extract all cells of the component and compute the feature
        
        components_cell_ids = vertcat(components_struct.cell.id);
        components_cell_surface_features = vertcat(components_struct.cell.surface_features);
        components_center = mean(components_cell_surface_features(:,1:3),1);
        [component_eigvector, ~, ~]  = eigenspace(components_cell_surface_features(:,1:3),0);

        % Extract boundary cell
        [components_bound_cell_ids, ~] = boundary_cells(Tree, components_cell_ids(:,1));

        % Searching cell connected to the boundary cells
        flag = false(numel(reference_cell_ids),1);
        for k = 1:numel(components_bound_cell_ids)
            neighbour_cell_ids = Window_Query_Neighbour_Cell(Tree, components_bound_cell_ids(k), reference_cell_ids, bridge_threshold.search_scale);
            mask = ismember(reference_cell_ids, neighbour_cell_ids);
            flag(mask) = true;
        end
        components_connect_bound_cell_ids = unique(reference_cell_ids(flag));
        components_connect_bound_cell_ids = setdiff(components_connect_bound_cell_ids,components_cell_ids(:,1));

        % Extract the cell containing the points of the parapet or safe guard:
        % on the bound or above 
        connect_component_region = struct('cell',[]);
        count = 0;
        % For connect_bound_cells
        for k = 1:numel(components_connect_bound_cell_ids)

            % Retrieve elevation from the component
            components_connect_bound_cell_id = components_connect_bound_cell_ids(k);
            neighbour_cell_ids = Window_Query_Neighbour_Cell(Tree, components_connect_bound_cell_id, components_cell_ids(:,1), bridge_threshold.search_scale);
            mask = ismember(components_cell_ids, neighbour_cell_ids);
            neighbour_cell_center_z = components_cell_surface_features(mask,3);

            % Retrieve the points in connect_bound-cell
            components_connect_bound_cell_ptc_ids = Tree.cell_pts(components_connect_bound_cell_id).id;
            components_connect_bound_cell_ptc = Tree.pts(components_connect_bound_cell_ptc_ids,1:3);

            % Extract points above the components
            mask = components_connect_bound_cell_ptc(:,3) >= mean(neighbour_cell_center_z) - 1.0*std(neighbour_cell_center_z);
            components_connect_bound_cell_ptc_ids = components_connect_bound_cell_ptc_ids(mask);

            if numel(components_connect_bound_cell_ptc_ids) >= threshold.min_num_pts
                components_connect_bound_cell_ptc = components_connect_bound_cell_ptc(mask,1:3);
                components_connect_bound_cell_height = max(components_connect_bound_cell_ptc(:,3)) - mean(neighbour_cell_center_z);

                % Assign to the connect_component_region
                connect_component_region.cell(count + 1).id = components_connect_bound_cell_id;
                connect_component_region.cell(count + 1).center = mean(components_connect_bound_cell_ptc,1);
                connect_component_region.cell(count + 1).ptc_ids = components_connect_bound_cell_ptc_ids;
                connect_component_region.cell(count + 1).height = components_connect_bound_cell_height;
                count = count + 1;
            end
        end

        % For component cell themselve
        for k = 1:size(components_cell_ids, 1)
            % Retrieve elevation from the component
            components_cell_id = components_cell_ids(k, 1);
            components_cell_ptc_ids = components_struct.cell(k).ptc_ids;
            components_cell_center_z = components_cell_surface_features(k,3);

            % Retrieve the points in connect_bound-cell
            components_cell_origin_ptc_ids = Tree.cell_pts(components_cell_id).id;
            components_cell_remain_ptc_ids = setdiff(components_cell_origin_ptc_ids,components_cell_ptc_ids);

            if numel(components_cell_remain_ptc_ids) >= threshold.min_num_pts
                % Extract origin points above the component
                components_cell_remain_ptc = Tree.pts(components_cell_remain_ptc_ids,1:3);
                mask = components_cell_remain_ptc(:,3) > components_cell_center_z - bridge_threshold.surface_error;
                components_cell_remain_ptc_ids = components_cell_remain_ptc_ids(mask);
                if numel(components_cell_remain_ptc_ids) >= threshold.min_num_pts
                    % Calculate the cell height
                    components_cell_remain_ptc = components_cell_remain_ptc(mask,1:3);
                    components_cell_height = max(components_cell_remain_ptc(:,3)) - components_cell_center_z;

                    % Assign to the connect_component_region
                    connect_component_region.cell(count + 1).id = components_cell_id;
                    connect_component_region.cell(count + 1).center = mean(components_cell_remain_ptc,1);
                    connect_component_region.cell(count + 1).ptc_ids = components_cell_remain_ptc_ids;
                    connect_component_region.cell(count + 1).height = components_cell_height;
                    count = count + 1;
                end
            end
        end

        % Use kernel estimation for extracting parapet and safe guard
        cell_height = vertcat(connect_component_region.cell.height);
        cell_center = vertcat(connect_component_region.cell.center);
        connect_component_features = zeros(numel(cell_height),2); %[id, height]
        [fi, xi, ~] = ksdensity(cell_height,'npoints', 100, 'bandwidth', bridge_threshold.bandwidth,'Kernel','epanechnikov');
        [~,locs] = findpeaks(fi, 'NPeaks', 3,'SortStr','descend');

    %     close all
    %     hold all
    %     plot(xi,fi,'b-')

        for k = 1:numel(locs)
            loc = locs(k);
            mask = (xi(loc) - bridge_threshold.bandwidth <= cell_height)&(cell_height <= xi(loc) + bridge_threshold.bandwidth);
            cell_height_peak = cell_height(mask);

    %         if bridge_threshold.parapet_height <= mean(cell_height_peak)
                connect_component_features(mask, 1) = k;
                connect_component_features(mask, 2) = cell_height_peak;
    %         end  
        end

        % Assign the cell size to the components
        component_region = struct('cell',[], 'status',[]);

        num_add_component = unique(connect_component_features(connect_component_features(:,1)>0, 1));
        for k = 1:size(num_add_component)
            local_ids = find(connect_component_features(:,1) == num_add_component(k));
            local_cell_center = cell_center(local_ids,1:3);
            proj_local_cell_center = proj_3Dpoints_3Dplane(local_cell_center, [components_center,component_eigvector(1,:)]);

            % Using clustering based on distance from the cell to center line of the
            % connection
            d_connection_component = dist_3Dpoints_3Dline(proj_local_cell_center, [components_center, component_eigvector(3,:)], component_eigvector(2,:));
    %         inlier_cell_ids = clustering_dbscan(d_connection_component, bridge_threshold.parapet_rail_cluster, bridge_threshold.min_no_cells, 1);
            inlier_cell_ids = clustering_dbscan(d_connection_component, threshold.cell_size, bridge_threshold.min_no_cells, 1);
            local_ids = local_ids(logical(inlier_cell_ids));

            component_region(k).status = 1;
            for count = 1:numel(local_ids)
                local_id = local_ids(count);
                component_region(k).cell(count).id = connect_component_region.cell(local_id).id;
                component_region(k).cell(count).ptc_ids = connect_component_region.cell(local_id).ptc_ids;
            end
        end

        % Update bridge data structure: child of child 
        bridge(bridge_tree_level).Link_List(:,end) = 0;
        for k = 1:length(component_region)
    %         bridge = bridge_tree_components(bridge, bridge_tree_level, parent_node, component_region, k, child_node + k);

            [bridge, ~, component_region] = bridge_tree_components(bridge, bridge_tree_level, new_parent_id, k, component_region);

        end

    end
end
%% Remove child
bridge(bridge_tree_level).Link_List(:,end) = [];
%% Write out the results
% ptc = [];
% for i = 1:length(Bridge(bridge_tree_level))
%     for j=1:length(Bridge(bridge_tree_level).Component)
%         ptc_ids = vertcat(Bridge(bridge_tree_level).Component(j).cell.ptc_ids);
%         ptc_xyz = OQTR.pts(ptc_ids,1:3);
%         ptc_xyz(:,4) = j;
%         ptc = [ptc;ptc_xyz];
%     end
%     
% end
% 
% our_dir = "c:\Users\ltruonghong\TUDelft\Point Cloud Processing\Matlab codes\Bridge segmentation\Test results_6\";
% outputFile = strcat(our_dir,'test.txt');
% write_txt(outputFile, ptc)

