function global_comp_link_list = connected_2D_components(Tree, cells_region, region_features, source_ids, target_ids, struct_threshold, threshold)
% This function is to determine connection between the region, which can be
% parallel or orthogonal. This based on an assumption the cell represent to the surface
% Input:
%     reference_id            : A id of the reference region
%     target_ids              : ids of the connected reion to the reference
%     scale                   : Searching window size
%     region_features         : Features of the region
%     source_ids              : reference region ids
%     target_ids              : A list of the region to be checked if they connect to source_id or not
%  Output:
%     component_link_list      : A list of the connection [ref_id, connected_ids]
%  Demo:
% 
%     cells_region
%     region_features
%     source_ids = bridge_bottom_surf_id;%bridge_deck_segment_ids
%     target_ids = connect_segs_ids;%candidate_components_segments_ids


%% Extract the boundary of the reference surface
global_comp_link_list = [];

while ~isempty(source_ids)
    % Retrieve boundary cells of the reference region
    source_id = source_ids(1);
    source_cell_ids = vertcat(cells_region(source_id).cell.id);
    [source_bound_cells_ids, ~] = get_boundary_cells(Tree, source_cell_ids(:,1));
    source_cell_surf_features = vertcat(cells_region(source_id).cell.surface_features);
%     source_cell_ptc_ids = vertcat(region(source_id).cell.ptc_ids);
%     source_cell_ptc_xyz = Tree.pts(source_cell_ptc_ids,1:3);
    
    % Retrieve features of the surface
    mask = ismember(region_features.ids, source_id);
    source_surf_cent = region_features.center(mask,:);
    source_surf_normal = region_features.normal(mask,:);
    source_surf_tangent = region_features.tangent(mask,:);
    
    % Search connection region if the number of target id is larger than 1
    if numel(target_ids) > 0
        % Remove source_ids out of target
        if ~isempty(global_comp_link_list)
            target_ids = setdiff(target_ids, reshape(global_comp_link_list(:,[1,2]),[],1),'stable');
        end
        
        % Remove source id if it available
        target_ids = setdiff(target_ids, source_id,'stable');
        
        % Preallocation
        comp_link_list = [inf,inf,inf];
        count = 1;
        for i = 1:numel(target_ids)

            % Extract cells on the bound
            cur_target_id = target_ids(i);
            cur_target_cell_ids = vertcat(cells_region(cur_target_id).cell.id);
            cur_target_cell_surf_feature = vertcat(cells_region(cur_target_id).cell.surface_features);
            [cur_target_edges_cell_ids, ~] = get_boundary_cells(Tree, cur_target_cell_ids(:,1));
        
            % Retrieve features of the surface
            mask = ismember(region_features.ids, cur_target_id);
            cur_target_surface_center = region_features.center(mask,:);
            cur_target_surface_normal = region_features.normal(mask,:);
    
%             cur_target_cell_ptc_ids = vertcat(region(cur_target_id).cell.ptc_ids);
%             current_target_cell_ptc_xyz = Tree.pts(cur_target_cell_ptc_ids,1:3);
%             
%             current_target_surface_center = mean(current_target_cell_ptc_xyz,1);
%             current_target_surface_normal = eigenspace(current_target_cell_ptc_xyz,1);
            
            % Check an intersection between two surface
            [intersect_line, check] = plane_plane_intersection([source_surf_cent, source_surf_normal], [cur_target_surface_center, cur_target_surface_normal], threshold.max_angle);
            
            % Establish connection
            check_connected_comp = false;
            if check == 1 % Intersection between two regions:
                
                % Determine the cells on region(s) close to the intersection line and projection of them
                % For a source region
                mask = ismember(source_cell_ids(:,1),source_bound_cells_ids); 
                source_bound_cells_cent = source_cell_surf_features(mask,1:3);
                dist_source_bound_cells_cent_intersect_line = dist_3Dpoints_3Dline(source_bound_cells_cent, intersect_line);
                mask = abs(dist_source_bound_cells_cent_intersect_line) <= threshold.cell_size;
                source_bound_cells_cent = source_bound_cells_cent(mask,:);
                
                % For target region
                mask = ismember(cur_target_cell_ids(:,1),cur_target_edges_cell_ids); 
                cur_target_bound_cells_cent = cur_target_cell_surf_feature(mask,1:3);
                dist_target_bound_cells_cent_intersect_line = dist_3Dpoints_3Dline(cur_target_bound_cells_cent, intersect_line);
                mask = abs(dist_target_bound_cells_cent_intersect_line) <= threshold.cell_size;
                cur_target_bound_cells_cent = cur_target_bound_cells_cent(mask,:);
                
                if (size(source_bound_cells_cent,1) > 5)&&(size(cur_target_bound_cells_cent,1) > 5)
                
                    % Find end points of source and target edges
                    [~, ~, source_bound_cells_end_cent, ~, ~] = findEndPts(source_bound_cells_cent);
                    source_edge_length = norm(source_bound_cells_end_cent(1:3) - source_bound_cells_end_cent(4:6));
                    [~, ~, cur_target_bound_cells_end_cent, ~, ~] = findEndPts(cur_target_bound_cells_cent);
                    cur_target_edge_length = norm(cur_target_bound_cells_end_cent(1:3) - cur_target_bound_cells_end_cent(4:6));
                    % Compute an overlap line
                    source_target_overlap_length = cal_line_segments_overlap(source_bound_cells_end_cent,cur_target_bound_cells_end_cent);
                    
                    if (source_target_overlap_length >= 0.5*source_edge_length) || (source_target_overlap_length >= 0.5*cur_target_edge_length)
                        check_connected_comp = true;  
                        
                        % Filter the region :
                    end
                end

            else % Parallel or adjoined surface
                % Find the connection between two regions: parallel case
                flag = false(numel(source_bound_cells_ids),1);
                for j = 1:numel(cur_target_edges_cell_ids)
                    neighbour_cell_ids = Window_Query_Neighbour_Cell(Tree, cur_target_edges_cell_ids(j), source_bound_cells_ids, struct_threshold.search_scale);
                    mask = (0 <= cur_target_cell_surf_feature(j,3) - source_cell_surf_features(ismember(source_cell_ids, neighbour_cell_ids),3))&...
                           (cur_target_cell_surf_feature(j,3) - source_cell_surf_features(ismember(source_cell_ids, neighbour_cell_ids),3) <= struct_threshold.curb_max_height);
                    neighbour_cell_ids = neighbour_cell_ids(mask);
                    mask = ismember(source_bound_cells_ids,neighbour_cell_ids); 
                    flag(mask) = true;
                end
                connected_cells_ids = source_bound_cells_ids(flag);
                clear j neighbour_cell_ids
                
                source_target_overlap_length = numel(connected_cells_ids)*threshold.cell_size;
                if source_target_overlap_length >= struct_threshold.adjoin_components_min_edge_length
                    check_connected_comp = true;
                end
            end
      
            % Check length of the connected cell
            if check_connected_comp
                comp_link_list(count,:) = [source_id, cur_target_id, check];
                count = count + 1;
            end
        end
        % Remove all preallocation link list
        mask = isinf(comp_link_list(:,1));
        comp_link_list(mask,:) = [];
    
        % Define a sign of the connection based on a sign distance & update a target region id
        if size(comp_link_list, 1) >= 1
            mask = ismember(region_features.ids, comp_link_list(:,2));
            comp_link_region_cent = region_features.center(mask,:);

            %  Divide into leave and right side
            comp_link_list(:,4) = ptc_sign_distance(comp_link_region_cent, [source_surf_cent, source_surf_tangent]);

            % Update a global link list
            global_comp_link_list = [global_comp_link_list; comp_link_list];

            clear component_link_region_center source_region_center source_tangent_vector
        end

        % Update source region id
        source_ids = union(source_ids, comp_link_list(:,2), 'stable');
    end
    
    % Update source_ids
    source_ids = setdiff(source_ids, source_id, 'stable');
    
end
%% 
  