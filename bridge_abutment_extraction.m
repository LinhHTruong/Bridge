function bridge_abutment = bridge_abutment_extraction(Tree, subTree, road_surface, struct_threshold, threshold, time_record)
%% This fuction is to segment the point cloud of the abutment
% Input:
%       Tree                    : Data structure stored entire data
%       subTree                 : Data structure stored 2D cells occupying the points of substructure
%       road_surface            : Data structure stored features of the road surface
%       struct_threshold
%       threshold         
% Output:
%       bridge_abutment         : Data structure stored points' indices of each plane

% Demo:
% Tree = OQTR;
% subTree = subOQTR;
% struct_threshold;
% threshold = THRESHOLD;
% road_surface = Road_Surface;
if time_record
    tic
end
%% Set up data structure
bridge_abutment = struct('plane', [], 'Description', {});
%% Retrieve cluster of the abutment
substruct_comp_name = vertcat(subTree.Description);
mask = cellfun(@(s) contains(s, 'Abutment', 'IgnoreCase', true), substruct_comp_name, 'UniformOutput', false);
abutment_ids = find(cell2mat(mask));

%% Set up the parameters for segmentation
voxel_region_growing_threshold = voxel_region_growing_params(1, struct_threshold.section_width, threshold.max_angle, threshold.sampling_step);

%% Voxel based segmentation
for i = 1:numel(abutment_ids)
%     i = 1;
    abutment_id = abutment_ids(i);
    % Retrieve the candidate points of the abutment
    comp_pts_ids = vertcat(subTree(abutment_id).cell.ptc_ids);
    comp_pts_xyz = Tree.pts(comp_pts_ids,1:3);
    
    % Voxel-based segmentation
    [ptc_segment_info, ~, ~] = voxel_region_growing_segmentation(comp_pts_xyz,voxel_region_growing_threshold);
    ptc_segment_info(:,1) = comp_pts_ids(ptc_segment_info(:,1)); % Convert to global indices
    
    % Component connectivity
    region_features = cal_region_features(Tree, ptc_segment_info);
    
    % Determine a stem: vertical and perpendicular to bridge deck
    road_surface_tangent = road_surface.tangent;
    road_surface_width = road_surface.mbb.short_edge;
    
    % Check based on an angle deviation and size
    region_proj_normal = vertcat(region_features.normal);
    region_proj_normal(:,3) = 0;
    road_surface_proj_tangent = road_surface_tangent;
    road_surface_proj_tangent(:,3) = 0;
    cosin_tangent_norm = abs(cosine_vectors(road_surface_proj_tangent, region_proj_normal));
    
    % Select first 5 segments having the small deviation: For the scase of
    % a skew bridge the stem may not be perpendicular to the road direction
    [~, maxk_ids] = maxk(cosin_tangent_norm,5);
    
    % Check dimensions of a stem
    mask = region_features.bounding_box(maxk_ids,2) >= 0.75*road_surface_width;
    maxk_ids = maxk_ids(mask);
    
    % Select the lowest segment
    [~, min_id] = min(region_features.center(maxk_ids,3));
    stem_surf_id = maxk_ids(min_id);
    
    % Establish a connectivity: The candidate segment must have a size
    % larger than the minimum size and their locations are not below the
    % stem
    mask = all((region_features.bounding_box(:,[2,3]) >= 0.5*struct_threshold.section_width),2)&...
               (region_features.bounding_box(:,1) >= 0.5*struct_threshold.section_width*struct_threshold.section_length)&...
           (region_features.center(:,3) >= region_features.center(stem_surf_id,3) - 0.5*region_features.bounding_box(stem_surf_id,3) + struct_threshold.section_width);
    connect_segs_ids = setdiff(region_features.ids(mask), stem_surf_id); 
    comp_link_list = connected_components_based_points(Tree, ptc_segment_info,...
                    region_features, stem_surf_id, connect_segs_ids, struct_threshold, threshold);

    % Filtering outlier points
    for j = 1:size(comp_link_list,1)
        ptc_segment_info = plane_outlier_removal(Tree, ptc_segment_info, comp_link_list(j,[1,2]), threshold);
    end

    % Set up the structure
    comp_region_ids = unique(reshape(comp_link_list(:,[1,2]),[],1));
    bridge_abutment(i).Description = strcat('Abutment' , " ", num2str(i));
    for j = 1:numel(comp_region_ids)
        mask = ptc_segment_info(:,2) == comp_region_ids(j);
        bridge_abutment(i).plane(j).ptc_ids =  ptc_segment_info(mask,1);
        
    end
end
if time_record
    fprintf('Executing time for extracting planes of the abutment: %.2f seconds \n', toc);
end
% % Write results of segmentation
% temp_ptc = [ptc(Segment_Ptc_Info(:,1),1:3), Segment_Ptc_Info(:,2)];
% seg_center = mean(ptc(Segment_Ptc_Info(:,1),1:3));
% outputFile = strcat('Abutment_', num2str(seg_center(1)),'_', num2str(seg_center(2)),'.txt');
% write_txt(outputFile, temp_ptc)    
% clear temp_ptc seg_center


% 
% %% Determine a stem: vertical and perpendicular to bridge deck
% bridge_deck_ptc_ids = vertcat(bridge(1).Component(1).cell.ptc_ids);
% bridge_deck_ptc_xyz = Tree.pts(bridge_deck_ptc_ids,1:3);
% [bridge_deck_eigen_vectors, ~, ~] = eigenspace(bridge_deck_ptc_xyz, 0);
% % [~,~,~,~,bridge_deck_dimensions] = min3Dboundingbox(bridge_deck_ptc_xyz(:,1),bridge_deck_ptc_xyz(:,2),bridge_deck_ptc_xyz(:,3),'v',3);
%     
% 
% % Condition 1: perpendicular to the vertical
% cosin_norm_oz = abs(cosine_vectors(segment_features(:,5:7), [0., 0., 1.]));
% 
% mask = cosin_norm_oz <= sin(pi/2 - deg2rad(2.0*threshold.max_angle));
% segment_features_condition_01 = segment_features(mask,:);
% 
% % Condition 2: parallel to the bridge deck
% cosin_norm_bridge_deck = abs(cosine_vectors(segment_features_condition_01(:,5:7), bridge_deck_eigen_vectors(3,:)));
% mask = cosin_norm_bridge_deck >= cos(deg2rad(threshold.max_angle));
% segment_features_condition_02 = segment_features_condition_01(mask,:);
% 
% % The bridge abutment stem is the largest 
% [~, max_id] = max(segment_features_condition_02(:,end));
% bridge_abutment_stem = segment_features_condition_02(max_id,1);
% 
% clear mask max_id bridge_deck_ptc_ids bridge_deck_ptc_xyz bridge_deck_eigen_vectors cosin_norm_oz cosin_norm_bridge_deck
% %% Searching the connection region
% % Check connectivity of the remain region -> need to be improved
% tic
% mask = ((segment_features(:,9) >= 0.2)|(segment_features(:,10) >= 0.2))&(segment_features(:,end) > 0.04);
% check_segment_ids = segment_features(mask,1);
% seeding_segment_ids = bridge_abutment_stem;
% % Preallocation component_link_list
% component_link_list = [0, seeding_segment_ids]; % Root
% 
% % Update check_segment
% mask = ismember(check_segment_ids, seeding_segment_ids);
% check_segment_ids(mask) = [];
% 
% % Searching merging region
% while ~isempty(seeding_segment_ids) 
% 
%     % Extract features of the seeding_region
%     current_seed_segment_id = seeding_segment_ids(1);
%     mask = ismember(segment_features(:,1),current_seed_segment_id);
%     current_seed_segment_features = segment_features(mask,:);
% 
%     % Retrieve features of check segment is
%     mask = ismember(segment_features(:,1),check_segment_ids);
%     check_segments_features = segment_features(mask,:);
%     
%     % Filtering -> select nearly perpendicular
%     seed_check_segment_cosin =  abs(cosine_vectors(current_seed_segment_features(:,5:7), check_segments_features(:,5:7)));
%     mask = seed_check_segment_cosin <= cos(pi/2 - deg2rad(2.0*threshold.max_angle));
%     check_segments_features = check_segments_features(mask,:);
%     
%     % Determine an intersection
%     flag = false(size(check_segments_features,1),1);
%     for j = 1:size(check_segments_features,1)
%         % Check if a pair of the segment is intersection
%         [intersection_line, check] = plane_plane_intersection(current_seed_segment_features(2:7), check_segments_features(j,2:7), threshold.max_angle);
% 
%         if check == 1
%             % Determine the voxels on segment(s) close to the intersection line and projection of them
%             % For a source region
%             current_seed_segment_voxel_ids = vertcat(Component_Region(current_seed_segment_id).voxel.id);
%             current_seed_segment_voxel_center = 0.5*(COCT.voxel_bounds(current_seed_segment_voxel_ids,4:6) + COCT.voxel_bounds(current_seed_segment_voxel_ids,1:3));
%             dist_seed_voxel_center_intersection_line = dist_3Dpoints_3Dline(current_seed_segment_voxel_center, intersection_line);
%             mask = abs(dist_seed_voxel_center_intersection_line) <= threshold.voxel_size;
%             current_seed_edges_voxel_center = current_seed_segment_voxel_center(mask,:);
% 
%             % For target region
%             current_check_segment_voxel_ids = vertcat(Component_Region(check_segments_features(j,1)).voxel.id);
%             current_check_segment_voxel_center = 0.5*(COCT.voxel_bounds(current_check_segment_voxel_ids,4:6) + COCT.voxel_bounds(current_check_segment_voxel_ids,1:3));
%             dist_check_voxel_center_intersection_line = dist_3Dpoints_3Dline(current_check_segment_voxel_center, intersection_line);
%             mask = abs(dist_check_voxel_center_intersection_line) <= threshold.voxel_size;
%             current_check_edges_voxel_center = current_check_segment_voxel_center(mask,:);
% 
%             if (size(current_seed_edges_voxel_center,1) > 1)&&(size(current_check_edges_voxel_center,1) > 15)
% 
%                 % Find end points of source and target edges
%                 [~, ~, current_seed_end_edges_voxel_center, ~, ~] = findEndPts(current_seed_edges_voxel_center);
%                 [~, ~, current_check_end_edges_voxel_center, ~, ~] = findEndPts(current_check_edges_voxel_center);
%                 current_seed_edge_length = norm(current_seed_end_edges_voxel_center(1:3) - current_seed_end_edges_voxel_center(4:6));
%                 current_check_edge_length = norm(current_check_end_edges_voxel_center(1:3) - current_check_end_edges_voxel_center(4:6));
%                 % Compute an overlap length of the line segme
%                 current_seed_check_overlap_length = cal_line_segments_overlap(current_seed_end_edges_voxel_center,current_check_end_edges_voxel_center);
%                 
%                 % Decision whether a real connection or not
%                 if ((current_seed_check_overlap_length >= 0.5*current_seed_edge_length) || (current_seed_check_overlap_length >= 0.5*current_check_edge_length))
%                      
%                     flag(j) = true;
%                 end
%                 clear current_seed_end_edges_voxel_center current_check_end_edges_voxel_center current_seed_edge_length current_check_edge_length current_seed_check_overlap_length  
%             end
%             clear mask current_seed_segment_voxel_ids current_seed_segment_voxel_center dist_seed_voxel_center_intersection_line current_seed_edges_voxel_center current_check_segment_voxel_center dist_check_voxel_center_intersection_line current_check_edges_voxel_center
%         end
%         
%     end
% 
%     % Update group region
%     connection_segment_ids = check_segments_features(flag, 1);
%     if ~isempty(connection_segment_ids)
%         % Update seeding region
%         mask = ismember(connection_segment_ids, component_link_list(:));
%         seeding_segment_ids = union(seeding_segment_ids, connection_segment_ids(mask));
%                 
%         % Connection list
%         [m,n] = ndgrid(current_seed_segment_id,connection_segment_ids);
%         new_component_link_list = [m(:),n(:)];
%         component_link_list = union(component_link_list, new_component_link_list, 'rows', 'stable');
%         
%     end
%     % remove the seed region to be check
%     seeding_segment_ids = setdiff(seeding_segment_ids, current_seed_segment_id);
% end  
% 
% fprintf('Running time for group region connectivity: %.2f seconds \n', toc);
% clear seeding_segment_ids current_seed_segment_id
% 
% %% Update bridge structure
% component_segments = struct('cell',[], 'status',[]);
% for i = 1:size(component_link_list,1)
%     component_id = component_link_list(i,2);
%     component_segments(i).status = 1;
%     for j = 1:length(Component_Region(component_id).voxel)
%         component_segments(i).cell(j).id = Component_Region(component_id).voxel(j).id; 
%         component_segments(i).cell(j).ptc_ids = ptc_ids(Component_Region(component_id).voxel(j).ptc_ids); % Transform from local to global
%     end
% end
%   
% for i = 1:length(component_segments)
% [bridge, ~, component_segments] = bridge_tree_components(bridge, bridge_tree_level, component_link_list(i,1), component_link_list(i,2), component_segments);
% end
% % remove old ids of the child components
% bridge(bridge_tree_level).Link_List(:,end) = [];

