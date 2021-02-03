function global_comp_link_list = connected_components_based_points(Tree, ptc_segment_info, region_features, source_ids, target_ids, struct_threshold, threshold)
% This function is to determine connection between the region, which can be
% parallel or orthogonal. This based on an assumption the cell represent to the surface
% Input:
%     reference_id            : A id of the reference region
%     target_ids              : ids of the connected reion to the reference
%     scale                   : Searching window size
%     Region_Features         : Features of the region
%     option:                 : parallel or orthorgonal
%     direction               : searching region connected to the reference region: long edge, short edge or all
%  Output:
%     component_link_list      : A list of the connection [ref_id, connected_ids]
%  Demo:
%     ptc_segment_info = ptc_segment_info;
%     region_features
%     source_ids = stem_surf_id;%bridge_deck_segment_ids
%     target_ids = connect_segs_ids;

%% Extract the boundary of the reference surface
global_comp_link_list = [];

while ~isempty(source_ids)
    % Retrieve boundary cells of the reference region
    source_id = source_ids(1);

    % Retrieve features of the surface
    mask = ismember(region_features.ids, source_id);
    source_surf_cent = region_features.center(mask,:);
    source_surf_normal = region_features.normal(mask,:);
%     source_surf_tangent = region_features.tangent(mask,:);
    
    % Remove source id if it available
    target_ids = setdiff(target_ids, source_id,'stable');
    
    % Remove segment asigned to the connected link out of target ids
    if ~isempty(global_comp_link_list)
        target_ids = setdiff(target_ids, reshape(global_comp_link_list(:,[1,2]),[],1),'stable');
    end

    % Search connection region if the number of target id is larger than 1
    if numel(target_ids) > 0
        
        % Preallocation
        comp_link_list = [inf,inf,inf];
        count = 1;
        for i = 1:numel(target_ids)

            % Retrieve features of the surface
            cur_target_id = target_ids(i);
            mask = ismember(region_features.ids, cur_target_id);
            target_surface_center = region_features.center(mask,:);
            target_surface_normal = region_features.normal(mask,:);
               
            % Check an intersection between two surface
            [intersect_line, check] = plane_plane_intersection([source_surf_cent, source_surf_normal], [target_surface_center, target_surface_normal], threshold.max_angle);
            
            % Establish connection
            check_connected_comp = false;
            if check == 1 % Intersection between two regions:
                
                % Determine edge length of each segment
                % For a source region
                mask = ismember(ptc_segment_info(:,2),source_id); 
                source_ptc_ids = ptc_segment_info(mask,1);
                source_ptc_xyz = Tree.pts(source_ptc_ids,1:3);
                
                % Project points on to its fitting surface
                source_proj_plane_ptc_xyz = proj_3Dptc_3Dplane(source_ptc_xyz,[source_surf_cent, source_surf_normal]);
                
                % Select the points within a buffer of the intersection line
                dist_source_line = dist_3Dpoints_3Dline(source_proj_plane_ptc_xyz, intersect_line);
                mask = abs(dist_source_line) <= struct_threshold.section_width;
                source_proj_plane_ptc_xyz = source_proj_plane_ptc_xyz(mask,:);
                
                % Project points onto the intersection line
                source_proj_line_ptc_xyz = proj_point_on_3Dline(source_proj_plane_ptc_xyz, intersect_line);
                
                % Find end points of source and target edges
                if isempty(source_proj_line_ptc_xyz)
                    continue;
                end
                [~, ~, source_edge_end_pts, ~, ~] = findEndPts(source_proj_line_ptc_xyz);
                source_edge_length = norm(source_edge_end_pts(1:3) - source_edge_end_pts(4:6));
                
                % For a target region
                mask = ismember(ptc_segment_info(:,2),cur_target_id); 
                target_ptc_ids = ptc_segment_info(mask,1);
                target_ptc_xyz = Tree.pts(target_ptc_ids,1:3);
                
                % Project points on to its fitting surface
                target_proj_plane_ptc_xyz = proj_3Dptc_3Dplane(target_ptc_xyz,[target_surface_center, target_surface_normal]);
                
                % Select the points within a buffer of the intersection line
                dist_target_line = dist_3Dpoints_3Dline(target_proj_plane_ptc_xyz, intersect_line);
                mask = abs(dist_target_line) <= struct_threshold.section_width;
                target_proj_plane_ptc_xyz = target_proj_plane_ptc_xyz(mask,:);
                
                % Project points onto the intersection line
                if isempty(target_proj_plane_ptc_xyz)
                    continue;
                end
                target_proj_line_ptc_xyz = proj_point_on_3Dline(target_proj_plane_ptc_xyz, intersect_line);
                
                % Find end points of source and target edges
                [~, ~, target_edge_end_pts, ~, ~] = findEndPts(target_proj_line_ptc_xyz);
                target_edge_length = norm(target_edge_end_pts(1:3) - target_edge_end_pts(4:6));
                  
                % Compute an overlap line
                source_target_overlap_length = cal_line_segments_overlap(source_edge_end_pts,target_edge_end_pts);
                % Chec if they intersection
                if (source_target_overlap_length >= 0.5*source_edge_length) || (source_target_overlap_length >= 0.5*target_edge_length)
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
    
        % Update a target region id
        if size(comp_link_list, 1) >= 1
            % Update a global link list
            global_comp_link_list = [global_comp_link_list; comp_link_list];
        end

        % Update source region id
        source_ids = union(source_ids, comp_link_list(:,2), 'stable');
    end
    
    % Update source_ids
    source_ids = setdiff(source_ids, source_id, 'stable');
    
end
%% 
  