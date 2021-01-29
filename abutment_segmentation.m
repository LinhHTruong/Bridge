function bridge = abutment_segmentation(Tree, ptc_ids, bridge, bridge_tree_level, threshold)

%% This fuction is to segment the point cloud of the abutment


% Demo:
% bridge = Bridge;
% bridge_tree_level = 3;
% Tree = OQTR;
% ptc_ids = substructure_cell_ptc_ids;
% threshold = THRESHOLD;
% bridge_threshold;
%% Retrieve data
ptc = Tree.pts(ptc_ids,1:3);
%% Voxel based segmentation
[Segment_Ptc_Info,Component_Region,COCT] = voxel_region_growing(ptc, threshold);%,'yes', out_file);

% % Write results of segmentation
% temp_ptc = [ptc(Segment_Ptc_Info(:,1),1:3), Segment_Ptc_Info(:,2)];
% seg_center = mean(ptc(Segment_Ptc_Info(:,1),1:3));
% outputFile = strcat('Abutment_', num2str(seg_center(1)),'_', num2str(seg_center(2)),'.txt');
% write_txt(outputFile, temp_ptc)    
% clear temp_ptc seg_center


%% Compute features for each segments
num_segment = max(Segment_Ptc_Info(:,2));
segment_features = zeros(num_segment, 11);

for i = 1:num_segment

    % Retrieve points of the segment
    mask = Segment_Ptc_Info(:,2) == i;
    segment_ptc_ids = Segment_Ptc_Info(mask,1);
    segment_ptc_xyz = ptc(segment_ptc_ids, 1:3);
    segment_center = mean(segment_ptc_xyz, 1);
    [segment_eigenvector, ~, segment_residual] = eigenspace(segment_ptc_xyz, 4);
   
    
%     
%     [~,~,~,~,segment_edgelength] = min3Dboundingbox(segment_ptc_xyz(:,1),segment_ptc_xyz(:,2),segment_ptc_xyz(:,3),'v',3);
%     if ~isempty(segment_edgelength)
%         segment_edgelength = sort(segment_edgelength, 'descend');
%         segment_surface_area = segment_edgelength(1)*segment_edgelength(2);
%     else
%         segment_edgelength = [-inf, -inf];
%         segment_surface_area = -inf;
%     end
%     segment_features(i,:) = [i,segment_center, segment_eigenvector(1,:),segment_residual, segment_edgelength(1:2)', segment_surface_area];    
   
    rot_euler_angle = vrrotvec(segment_eigenvector(1,:),[0.0,0.0,1.0]);
    rot_matrix = vrrotvec2mat(rot_euler_angle);
    
    rot_segment_ptc_xyz = segment_ptc_xyz*rot_matrix;
    segment_min_bb = min2DBoundingBox(rot_segment_ptc_xyz(:,1:2)');
    segment_features(i,:) = [i,segment_center, segment_eigenvector(1,:),segment_residual, segment_min_bb.long_edge, segment_min_bb.short_edge, segment_min_bb.area];    
   
    
     clear mask segment_ptc_ids segment_ptc_xyz segment_center segment_eigenvector segment_edgelength segment_surface_area
end
%% Determine a stem: vertical and perpendicular to bridge deck
bridge_deck_ptc_ids = vertcat(bridge(1).Component(1).cell.ptc_ids);
bridge_deck_ptc_xyz = Tree.pts(bridge_deck_ptc_ids,1:3);
[bridge_deck_eigen_vectors, ~, ~] = eigenspace(bridge_deck_ptc_xyz, 0);
% [~,~,~,~,bridge_deck_dimensions] = min3Dboundingbox(bridge_deck_ptc_xyz(:,1),bridge_deck_ptc_xyz(:,2),bridge_deck_ptc_xyz(:,3),'v',3);
    

% Condition 1: perpendicular to the vertical
cosin_norm_oz = abs(cosine_vectors(segment_features(:,5:7), [0., 0., 1.]));

mask = cosin_norm_oz <= sin(pi/2 - deg2rad(2.0*threshold.max_angle));
segment_features_condition_01 = segment_features(mask,:);

% Condition 2: parallel to the bridge deck
cosin_norm_bridge_deck = abs(cosine_vectors(segment_features_condition_01(:,5:7), bridge_deck_eigen_vectors(3,:)));
mask = cosin_norm_bridge_deck >= cos(deg2rad(threshold.max_angle));
segment_features_condition_02 = segment_features_condition_01(mask,:);

% The bridge abutment stem is the largest 
[~, max_id] = max(segment_features_condition_02(:,end));
bridge_abutment_stem = segment_features_condition_02(max_id,1);

clear mask max_id bridge_deck_ptc_ids bridge_deck_ptc_xyz bridge_deck_eigen_vectors cosin_norm_oz cosin_norm_bridge_deck
%% Searching the connection region
% Check connectivity of the remain region -> need to be improved
tic
mask = ((segment_features(:,9) >= 0.2)|(segment_features(:,10) >= 0.2))&(segment_features(:,end) > 0.04);
check_segment_ids = segment_features(mask,1);
seeding_segment_ids = bridge_abutment_stem;
% Preallocation component_link_list
component_link_list = [0, seeding_segment_ids]; % Root

% Update check_segment
mask = ismember(check_segment_ids, seeding_segment_ids);
check_segment_ids(mask) = [];

% Searching merging region
while ~isempty(seeding_segment_ids) 

    % Extract features of the seeding_region
    current_seed_segment_id = seeding_segment_ids(1);
    mask = ismember(segment_features(:,1),current_seed_segment_id);
    current_seed_segment_features = segment_features(mask,:);

    % Retrieve features of check segment is
    mask = ismember(segment_features(:,1),check_segment_ids);
    check_segments_features = segment_features(mask,:);
    
    % Filtering -> select nearly perpendicular
    seed_check_segment_cosin =  abs(cosine_vectors(current_seed_segment_features(:,5:7), check_segments_features(:,5:7)));
    mask = seed_check_segment_cosin <= cos(pi/2 - deg2rad(2.0*threshold.max_angle));
    check_segments_features = check_segments_features(mask,:);
    
    % Determine an intersection
    flag = false(size(check_segments_features,1),1);
    for j = 1:size(check_segments_features,1)
        % Check if a pair of the segment is intersection
        [intersection_line, check] = plane_plane_intersection(current_seed_segment_features(2:7), check_segments_features(j,2:7), threshold.max_angle);

        if check == 1
            % Determine the voxels on segment(s) close to the intersection line and projection of them
            % For a source region
            current_seed_segment_voxel_ids = vertcat(Component_Region(current_seed_segment_id).voxel.id);
            current_seed_segment_voxel_center = 0.5*(COCT.voxel_bounds(current_seed_segment_voxel_ids,4:6) + COCT.voxel_bounds(current_seed_segment_voxel_ids,1:3));
            dist_seed_voxel_center_intersection_line = dist_3Dpoints_3Dline(current_seed_segment_voxel_center, intersection_line);
            mask = abs(dist_seed_voxel_center_intersection_line) <= threshold.voxel_size;
            current_seed_edges_voxel_center = current_seed_segment_voxel_center(mask,:);

            % For target region
            current_check_segment_voxel_ids = vertcat(Component_Region(check_segments_features(j,1)).voxel.id);
            current_check_segment_voxel_center = 0.5*(COCT.voxel_bounds(current_check_segment_voxel_ids,4:6) + COCT.voxel_bounds(current_check_segment_voxel_ids,1:3));
            dist_check_voxel_center_intersection_line = dist_3Dpoints_3Dline(current_check_segment_voxel_center, intersection_line);
            mask = abs(dist_check_voxel_center_intersection_line) <= threshold.voxel_size;
            current_check_edges_voxel_center = current_check_segment_voxel_center(mask,:);

            if (size(current_seed_edges_voxel_center,1) > 1)&&(size(current_check_edges_voxel_center,1) > 15)

                % Find end points of source and target edges
                [~, ~, current_seed_end_edges_voxel_center, ~, ~] = findEndPts(current_seed_edges_voxel_center);
                [~, ~, current_check_end_edges_voxel_center, ~, ~] = findEndPts(current_check_edges_voxel_center);
                current_seed_edge_length = norm(current_seed_end_edges_voxel_center(1:3) - current_seed_end_edges_voxel_center(4:6));
                current_check_edge_length = norm(current_check_end_edges_voxel_center(1:3) - current_check_end_edges_voxel_center(4:6));
                % Compute an overlap length of the line segme
                current_seed_check_overlap_length = cal_line_segments_overlap(current_seed_end_edges_voxel_center,current_check_end_edges_voxel_center);
                
                % Decision whether a real connection or not
                if ((current_seed_check_overlap_length >= 0.5*current_seed_edge_length) || (current_seed_check_overlap_length >= 0.5*current_check_edge_length))
                     
                    flag(j) = true;
                end
                clear current_seed_end_edges_voxel_center current_check_end_edges_voxel_center current_seed_edge_length current_check_edge_length current_seed_check_overlap_length  
            end
            clear mask current_seed_segment_voxel_ids current_seed_segment_voxel_center dist_seed_voxel_center_intersection_line current_seed_edges_voxel_center current_check_segment_voxel_center dist_check_voxel_center_intersection_line current_check_edges_voxel_center
        end
        
    end

    % Update group region
    connection_segment_ids = check_segments_features(flag, 1);
    if ~isempty(connection_segment_ids)
        % Update seeding region
        mask = ismember(connection_segment_ids, component_link_list(:));
        seeding_segment_ids = union(seeding_segment_ids, connection_segment_ids(mask));
                
        % Connection list
        [m,n] = ndgrid(current_seed_segment_id,connection_segment_ids);
        new_component_link_list = [m(:),n(:)];
        component_link_list = union(component_link_list, new_component_link_list, 'rows', 'stable');
        
    end
    % remove the seed region to be check
    seeding_segment_ids = setdiff(seeding_segment_ids, current_seed_segment_id);
end  

fprintf('Running time for group region connectivity: %.2f seconds \n', toc);
clear seeding_segment_ids current_seed_segment_id

%% Update bridge structure
component_segments = struct('cell',[], 'status',[]);
for i = 1:size(component_link_list,1)
    component_id = component_link_list(i,2);
    component_segments(i).status = 1;
    for j = 1:length(Component_Region(component_id).voxel)
        component_segments(i).cell(j).id = Component_Region(component_id).voxel(j).id; 
        component_segments(i).cell(j).ptc_ids = ptc_ids(Component_Region(component_id).voxel(j).ptc_ids); % Transform from local to global
    end
end
  
for i = 1:length(component_segments)
[bridge, ~, component_segments] = bridge_tree_components(bridge, bridge_tree_level, component_link_list(i,1), component_link_list(i,2), component_segments);
end
% remove old ids of the child components
bridge(bridge_tree_level).Link_List(:,end) = [];

