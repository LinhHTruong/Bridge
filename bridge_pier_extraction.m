function bridge_pier = bridge_pier_extraction(Tree, subTree, struct_threshold, threshold, time_record)
%% This fuction is to segmentation the point cloud of the pier
% Input:
%       Tree                    : Data structure stored entire data
%       subTree                 : Data structure stored 2D cells occupying the points of substructure
%       struct_threshold
%       threshold         

% Output:
%       bridge_pier         : Data structure stored points' indices of each plane of the pier

% Demo:
% Tree = OQTR;
% subTree = subOQTR;
% struct_threshold;
% threshold = THRESHOLD;
% road_surface = Road_Surface;

%% Set up data structure
bridge_pier = struct('plane', [], 'Description', {});

%% Retrieve cluster of the pier
substruct_comp_name = vertcat(subTree.Description);
mask = cellfun(@(s) contains(s, 'Pier', 'IgnoreCase', true), substruct_comp_name, 'UniformOutput', false);
pier_ids = find(cell2mat(mask));

%% Check types of the pier
for i = 1:numel(pier_ids)
%     i = 1;
    pier_id = pier_ids(i);
    % Retrieve the candidate points of the abutment
    comp_ptc_ids = vertcat(subTree(pier_id).cell.ptc_ids);
    comp_ptc_xyz = Tree.pts(comp_ptc_ids,1:3);
    
    % Extract the component of the pier
    pier_clusters = bridge_pier_type(comp_ptc_xyz, struct_threshold);

    % Segmentation
    for j = 1:size(pier_clusters)
        % Retrieve the points of the segment
        s_elev = pier_clusters(j,1) - struct_threshold.section_interval;
        e_elev = pier_clusters(j,2) + struct_threshold.section_interval;
        mask = (s_elev <= comp_ptc_xyz(:,3))&(comp_ptc_xyz(:,3) <= e_elev);
        pier_cluster_ptc_ids = comp_ptc_ids(mask);
        pier_cluster_ptc_xyz = comp_ptc_xyz(mask,:);
       
        % Set up parameters for the segmentation
        voxel_region_growing_threshold = voxel_region_growing_params(0, 0.5*pier_clusters(j,3), 2.0*threshold.max_angle, threshold.sampling_step);
        % Segmentation
        tic
        [ptc_segment_info, ~, ~] = voxel_region_growing_segmentation(pier_cluster_ptc_xyz,voxel_region_growing_threshold);
        ptc_segment_info(:,1) = pier_cluster_ptc_ids(ptc_segment_info(:,1)); % Update a global indices of the points
        toc
        % Filtering the 
        if pier_clusters(j,end) == 1
            % This is a pier cap
            
        elseif pier_clusters(j,end) == 2
            % They are pile cap
        else
            % Column
            no_columns = pier_count_column(pier_cluster_ptc_xyz, struct_threshold);
            
            % Compute the region features
            segs_features = cal_3D_region_feature(Tree, ptc_segment_info);
            
            % Set up the link of segments parts of the column
            % Preallocation 
            segs_link_list = inf(numel(segs_features.ids),2); %seg_ids, col_id]
            segs_link_list(:,1) = segs_features.ids;
            % Initialize
            flag_seg = true(numel(segs_features.ids),1);
            
            col_id = 1;
            while any(flag_seg)
                % Establish the connection link
                ref_id = segs_features.ids(find(flag_seg, 1, 'first'));
                target_ids = setdiff(segs_features.ids(flag_seg), ref_id);
                if ~isempty(target_ids)
                    link_list = connected_components_based_points(Tree, ptc_segment_info,...
                                        segs_features, ref_id, target_ids, struct_threshold, threshold);
                else
                    link_list = [ref_id, ref_id];
                end
                % Assign the segs_linked_list
                link_list = unique(reshape(link_list(:,[1,2]),[],1));
                mask = ismember(segs_link_list(:,1), link_list); 
                segs_link_list(mask,2) = col_id;
                col_id = col_id + 1;
                
                % Inactive segments assigned to the link
                flag_seg(mask) = false;
            end
            
            % Assign the column to the pier
            for k = 1:no_columns
                % Retrieve the segments of the column
                mask = segs_link_list(:,2) == k;
                seg_ids = segs_link_list(mask,1);
                
                bridge_pier(i).Description = strcat('Column', " ", num2str(k));
                for count = 1:numel(seg_ids)
                    % Retrieve the points in the segments
                    mask = ismember(ptc_segment_info(:,2), seg_ids(count));
                    plane_ptc_ids = ptc_segment_info(mask,1);

                    % Assign to the pier
                    bridge_pier(i).plane(count).ptc_ids = plane_ptc_ids;
                end
            end
        end

    end
    
end  

% %% Down sampling
% % Compute the features of the segment
% region_features = cal_region_features(Tree, ptc_segment_info);
% 
% % Filter based on the angle
% 
% 
% 
% 
% 
% %% Segmentation -> Identify columns
% tic
% kdtree = KDTreeSearcher(Section.center, "BucketSize",3);
% knn_all_ptc_ids = knnsearch(kdtree, Section.center, "k",size(Section.center,1));
% global_region = zeros(numel(Section.id),2);
% global_region(:,1) = Section.id;
% check_region_ids = Section.id;
% region_no = 0;
% while numel(check_region_ids) > 1%~isempty(check_region_ids)
% 
% %     if numel(check_region_ids) == 1
% %         region_no = region_no + 1;
% %         mask = global_region(:,1) == check_region_ids;
% %         global_region(mask,2) = region_no;
% %         check_region_ids = [];
% %         break;
% %     end
%     % Find seeding region
%     [~, min_id] = min(Section.center(check_region_ids,3));
%     start_seeding_ptc_ids = check_region_ids(min_id);
%     
%     % Find another point of the line segment
%     knn_ptc_ids = knnsearch(kdtree,Section.center(start_seeding_ptc_ids,1:3),'k',2);
%     end_seeding_ptc_ids = setdiff(knn_ptc_ids, start_seeding_ptc_ids);
%     
%     % Compute a horizontal distance
%     dist_ptc_ptc = norm(Section.center(start_seeding_ptc_ids,1:2) - Section.center(end_seeding_ptc_ids,1:2));
%     if dist_ptc_ptc <= pier.segment_distance
%         % Create a current segment
%         current_region_ptc_ids = [start_seeding_ptc_ids, end_seeding_ptc_ids];
%         % Update check column_ptc ids
%         check_region_ids = setdiff(check_region_ids,current_region_ptc_ids);
%     else
%         % No segment to be esablish -> start_seeding_ptc_ids: define a new
%         % segment
%         region_no = region_no + 1;
%         mask = global_region(:,1) == start_seeding_ptc_ids;
%         global_region(mask,2) = region_no;
%         check_region_ids = setdiff(check_region_ids, start_seeding_ptc_ids);
%     end
%     
%     % Region growing process
%     while (~isempty(end_seeding_ptc_ids))&&(~isempty(check_region_ids))
%         % Compute a region features -> fit line
%         current_region_center = mean(Section.center(current_region_ptc_ids,:),1);
%         [current_region_tangent, ~, ~] = eigenspace(Section.center(current_region_ptc_ids,:),2);
%         
%         % Find the close points
%         neighbour_ptc_ids = knn_all_ptc_ids(end_seeding_ptc_ids,2:end);
%         mask = ismember(neighbour_ptc_ids,check_region_ids);
%         neighbour_ptc_ids = neighbour_ptc_ids(mask);
%         % Compute distance from the points to the region
%         dist_ptc_line = distancePointLine3d(Section.center(neighbour_ptc_ids,:), [current_region_center, current_region_tangent]);
%         [min_val, min_id] = min(dist_ptc_line);
%         if min_val <= pier.segment_distance
%             % Update the end points
%             end_seeding_ptc_ids = neighbour_ptc_ids(min_id);
%             % Update the current region 
%             current_region_ptc_ids = union(current_region_ptc_ids,end_seeding_ptc_ids);
%             
%             % Update check region
%             check_region_ids = setdiff(check_region_ids, end_seeding_ptc_ids);
%         else
%             end_seeding_ptc_ids = [];
%         end
%     end
%     
%     % Update global region
%     region_no = region_no + 1;
%     mask = ismember(global_region(:,1),current_region_ptc_ids);
%     global_region(mask, 2) = region_no;   
% end
% toc
% %% Filtering the column
% tic
% Pier = struct("Component",[], "ptc_ids", [], "status", []);
% count = length(Pier);
% if count ~= 1
%     count = count + 1;
% end
% num_region = unique(global_region(global_region(:,2)>0,2));
% for i = 1:numel(num_region)
%     % Retrieve sections within the segment and filter outlier
%     mask = global_region(:,2) == i;
%     column_region_section_ids = global_region(mask,1);
%     column_region_section_boundingbox = Section.bounding_box(column_region_section_ids,1:2);
%     db_ids = dbscan(column_region_section_boundingbox(:,1:2),pier.section_bbox_error,pier.min_no_section);
%     mask = db_ids >= 0;
%     db_ids = db_ids(mask);
%     column_region_section_ids = column_region_section_ids(mask);
%     % Find the largest segment
%     num_segment = max(db_ids);
%     segment_count = histcounts(db_ids, num_segment);
%     segment_stat(:, 1) = unique(db_ids);
%     segment_stat(:, 2) = segment_count;
%     [~, ids] = max(segment_stat(:, 2));
%     largest_segment_id = segment_stat(ids,:);
%     
%     % Update largest segment
%     mask = db_ids==largest_segment_id;
%     column_region_section_ids = column_region_section_ids(mask);
%     column_region_section_boundingbox = Section.bounding_box(column_region_section_ids,1:2);
%     
%     % Extract the points along the section
%     column_region_section_center = mean(Section.center(column_region_section_ids,:),1);
%     [column_region_section_tangent, ~, ~] = eigenspace(Section.center(column_region_section_ids,:),2);
%     dist_ptc_column_center = distancePointLine3d(pier_ptc_xyz, [column_region_section_center, column_region_section_tangent]);
%     mask = dist_ptc_column_center <= 0.5*max(column_region_section_boundingbox(:)) + pier.section_bbox_error;
%     
%     pier_column_candidate_ptc_ids = pier_ptc_ids(mask);
%     pier_column_candidate_ptc_xyz = pier_ptc_xyz(mask,:);
% %     column_candidate_ptc_global_ids = ptc_ids(column_candidate_ptc_local_ids);
%     % Segmentation
%     [Segment_Ptc_Info,Component_Region,CTree] = voxel_region_growing(pier_column_candidate_ptc_xyz,threshold);
%     
%     % Write results of segmentation
%     temp_ptc = [pier_column_candidate_ptc_xyz(Segment_Ptc_Info(:,1),1:3), Segment_Ptc_Info(:,2)];
%     seg_center = mean(pier_column_candidate_ptc_xyz(Segment_Ptc_Info(:,1),1:3));
%     outputFile = strcat('Column_', num2str(seg_center(1)),'_', num2str(seg_center(2)),'.txt');
%     write_txt(outputFile, temp_ptc)    
%     clear temp_ptc seg_center
% 
%     % Extract the segment of the columns
%     [pier_column_ptc_local_ids, ~] = merging_column_segment(Segment_Ptc_Info, Component_Region, CTree, threshold);
%     
%     % Update data structure
%     Pier(count).Component = strcat("Column",num2str(i));
%     Pier(count).status = 1;
%     Pier(count).ptc_ids = pier_column_candidate_ptc_ids(pier_column_ptc_local_ids);% ptc_global_ids;
%     count = count + 1;
% end
% toc
% 
% % Not use for Germany bridge
% % %% Extract components of the pier caps
% % pier_column_ptc_ids = vertcat(Pier.ptc_ids);
% % pier_cap_candidate_ptc_ids = setdiff(pier_ptc_ids, pier_column_ptc_ids);
% % pier_cap_candidate_ptc_xyz = Tree.pts(pier_cap_candidate_ptc_ids,1:3);
% % 
% % % Extract the candidate points for a pier cap
% % mask = pier_cap_candidate_ptc_xyz(:,3) >= mean(Tree.pts(pier_column_ptc_ids,3));
% % pier_cap_candidate_ptc_ids = pier_cap_candidate_ptc_ids(mask);
% % pier_cap_candidate_ptc_xyz = pier_cap_candidate_ptc_xyz(mask,1:3);
% % 
% % 
% % %% Use voxel-based for large surfaces
% % [~,Component_Region,~] = voxel_region_growing(pier_cap_candidate_ptc_xyz,threshold);
% % 
% % %% Filtering the small segment
% % tic
% % 
% % num_segment = length(Component_Region);
% % segments_features = zeros(num_segment, 11);
% % 
% % for i = 1:num_segment
% %     % Retrieve points of the segment
% %     segment_ptc_xyz = vertcat(Component_Region(i).voxel.ptc_xyz);
% %     segment_center = mean(segment_ptc_xyz, 1);
% %     [segment_eigenvector, ~, segment_residual] = eigenspace(segment_ptc_xyz, 4);
% %     
% %     [~,~,~,~,segment_edgelength] = min3Dboundingbox(segment_ptc_xyz(:,1),segment_ptc_xyz(:,2),segment_ptc_xyz(:,3),'v',3);
% %     if ~isempty(segment_edgelength)
% %         segment_edgelength = sort(segment_edgelength, 'descend');
% %         segment_surface_area = segment_edgelength(1)*segment_edgelength(2);
% %     else
% %         segment_edgelength = [-inf, -inf];
% %         segment_surface_area = -inf;
% %     end
% %     segments_features(i,:) = [i,segment_center, segment_eigenvector(1,:),segment_residual, segment_edgelength(1:2), segment_surface_area];    
% %     clear mask segment_ptc_ids segment_ptc_xyz segment_center segment_eigenvector segment_edgelength segment_surface_area
% % end
% % 
% % mask = (segments_features(:,end) >= pier.component_min_area)& all(~isinf(segments_features(:,2)), 2);
% % segments_features = segments_features(mask,:);
% % 
% % %% Set up a connectivity -> continue work with other bridge
% % if ~isempty(segments_features)
% %     
% %     
% %     Pier_Cap_Component = struct('id',[], 'ptc_ids',[],'ptc_xyz',[]);
% %     mask = (segments_features(:,3) >= pier.component_min_area)&(~isinf(segments_features(:,2)));
% %     pier_cap_component_ids = segments_features(mask,1);
% %     pier_cap_non_component_ids = segments_features(~mask,1);
% %     for i = 1:numel(pier_cap_component_ids)
% %         Pier_Cap_Component(i).id = i;
% %         Pier_Cap_Component(i).ptc_ids = pier_cap_candidate_ptc_ids(vertcat(Component_Region(pier_cap_component_ids(i)).voxel.ptc_ids));
% %         Pier_Cap_Component(i).ptc_xyz = vertcat(Component_Region(pier_cap_component_ids(i)).voxel.ptc_xyz);
% %     end
% % end
% % 
% % toc
% %% Update bridge structure
% component_segments = struct('cell',[], 'status',[]);
% for i = 1:length(Pier)
%     component_segments(i).status = 1;
%     component_segments(i).cell.id = i;
%     component_segments(i).cell.ptc_ids = Pier(i).ptc_ids;
% end
%   
% for i = 1:length(Pier)
% [bridge, ~, component_segments] = bridge_tree_components(bridge, bridge_tree_level, 0, i, component_segments);
% end
% bridge(bridge_tree_level).Link_List(:,end) = [];
% 
% 
% 
% 
% 
% % 
% % 
% % 
% % 
% % %% remove old ids of the child components
% % 
% % 
% % 
% % 
% % bridge = Bridge;
% % bridge_tree_level = length(Bridge) + 1;
% % 
% % component_segments = struct('cell',[], 'status',[]);
% % for i = 1:size(component_link_list,1)
% %     component_id = component_link_list(i,2);
% %     component_segments(i).status = 1;
% %     for j = 1:length(Component_Region(component_id).voxel)
% %         component_segments(i).cell(j).id = Component_Region(component_id).voxel(j).id; 
% %         component_segments(i).cell(j).ptc_ids = ptc_ids(Component_Region(component_id).voxel(j).ptc_ids); % Transform from local to global
% %     end
% % end
% % %% Extracting the small surfaces by using point-based region growing
% % % pier_cap_component_ptc_ids = (arrayfun(@(x) vertcat(Component_Region(x).voxel.ptc_ids),pier_cap_component_ids,'UniformOutput',false));
% % % pier_cap_component_ptc_ids = cell2mat(pier_cap_component_ptc_ids); 
% % 
% % % Remaining points
% % pier_cap_non_component_ptc_ids = (arrayfun(@(x) vertcat(Component_Region(x).voxel.ptc_ids),pier_cap_non_component_ids,'UniformOutput',false));
% % pier_cap_non_component_ptc_ids = cell2mat(pier_cap_non_component_ptc_ids); 
% %  
% % 
% % out_file = strcat(our_dir, strtok(file_name,'.'),'_Pier_03_', 'pier_cap_Segmentation_all2.txt');
% % 
% % Component_Point_Segment_Info = PointBasedRegionGrowing(pier_cap_candidate_ptc_xyz(pier_cap_non_component_ptc_ids,1:3),threshold_ptc_segment,'y', out_file);
% % 
% 
% %% merge segment of the pier-column
% function [pier_column_ptc_ids, pier_column_segment_ids] = merging_column_segment(segment_ptc_info, Component_Region, CTree, threshold)
%     % This function is to merging segments of the columns
% 
%     % Input:
%     % Output:
%     % Demo;
%     %     segment_ptc_info = Segment_Ptc_Info;
% 
%     %
%     % Compute the feature of the columns
%     num_segments = max(segment_ptc_info(:,end));
%     all_segment_features = zeros(num_segments,2);
%     for i = 1:num_segments
%         % Compute length of the segment
%         mask = segment_ptc_info(:,2) == i;
%         segment_ptc_ids = segment_ptc_info(mask,1);
%         segment_ptc_z = CTree.pts(segment_ptc_ids, 3);
%         segment_length = max(segment_ptc_z) - min(segment_ptc_z);
%         all_segment_features(i,:) = [i, segment_length]; 
% 
%         clear mask segment_ptc_ids segment_ptc_z segment_length 
%     end
% 
%     % Find the largest segment connected to similar segments
%     % Preallocation 
%     seeding_segment_ids = 1; % the largest segment
%     pier_column_segment_ids = seeding_segment_ids; % Root
% 
%     while ~isempty(seeding_segment_ids)
%         % Retrive data of the seed segment
%         current_seed_segment_id = seeding_segment_ids(1);
%         mask = all_segment_features(:,1) == current_seed_segment_id;
%         current_seed_segment_length = all_segment_features(mask, 2);
% 
%         % Retrive data of the check segments
%         mask = ismember(all_segment_features(:,1),pier_column_segment_ids);
%         check_segment_features = all_segment_features(~mask, :);
% 
% 
%         % Check length
%         mask = (0.75*current_seed_segment_length <= check_segment_features(:,2))&(check_segment_features(:,2) <= 1.25*current_seed_segment_length);
%         check_segment_features = check_segment_features(mask,:);
% 
%         if ~isempty(check_segment_features)
%             % Retrieve voxel id in the current seed segment
%             current_seed_segment_voxel_ids = vertcat(Component_Region(current_seed_segment_id).voxel.id);
% 
%             % Preallocation
%             flag = false(size(check_segment_features,1),1);
%             for i = 1:size(check_segment_features,1)
% 
%                 % Retrieve voxel id in the current check segment
%                 current_check_segment_voxel_ids = vertcat(Component_Region(check_segment_features(i,1)).voxel.id);
% 
%                 % Check intersection: two segment share voxel
%                 seed_check_intersect_voxel_ids = intersect(current_seed_segment_voxel_ids,current_check_segment_voxel_ids);
%                 if numel(seed_check_intersect_voxel_ids)*threshold.voxel_size >= 0.5*current_seed_segment_length
%                     flag(i) = true;
%                 else
%                     % Two segment adjacent
%                     seed_check_neighbour_voxel_ids = query26_neighbour_voxels(CTree, current_seed_segment_voxel_ids, current_check_segment_voxel_ids);
%                     if numel(seed_check_neighbour_voxel_ids)*threshold.voxel_size >= 0.5*current_seed_segment_length
%                         flag(i) = true;
%                     else
%                         check_seed_neighbour_voxel_ids = query26_neighbour_voxels(CTree, current_check_segment_voxel_ids, current_seed_segment_voxel_ids);
%                         if numel(check_seed_neighbour_voxel_ids)*threshold.voxel_size >= 0.5*current_seed_segment_length
%                             flag(i) = true;
%                         end
%                     end
%                 end
%             end
% 
%             add_segment_ids = check_segment_features(flag,1);
% 
%             if ~isempty(add_segment_ids)
%                 % Update connectivity link list
%                 pier_column_segment_ids = union(pier_column_segment_ids,add_segment_ids,'stable');
% 
%                 % Update seed segment
%                 seeding_segment_ids = union(seeding_segment_ids,add_segment_ids,'stable');
%             end
% 
%         end
% 
%         % Update seed segment
%         seeding_segment_ids = setdiff(seeding_segment_ids, current_seed_segment_id, 'stable');
%     end
% 
%     % Return results
%     mask = ismember(segment_ptc_info(:,2), pier_column_segment_ids);
%     pier_column_ptc_ids = segment_ptc_info(mask, 1);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% % outputFile = strcat(our_dir, strtok(filename,'.'),'_Pier_test.txt');
% % write_txt(outputFile, section_ptc_xyz)
% % 
% % outputFile = strcat(our_dir, strtok(filename,'.'),'_Pier_Section_pos.txt');
% % write_txt(outputFile, [Section.center, global_region(:,2)])
% % 
% % outputFile = strcat(our_dir, strtok(filename,'.'),'_Pier_test_column_1.txt');
% % write_txt(outputFile, column_ptc_xyz)
% % 
% % outputFile = strcat(our_dir, strtok(filename,'.'),'_Pier_test_remain_ptc.txt');
% % write_txt(outputFile, remain_ptc_xyz)
% 
% 
% 
% 
% %% 
% 
% % threshold_ptc_segment.angle = 5;
% % threshold_ptc_segment.residual = 0.005;
% % threshold_ptc_segment.distance = 0.005;
% % threshold_ptc_segment.min_ptc = 5;
% % threshold_ptc_segment.cosin_angle = cos(threshold_ptc_segment.angle*pi/180); %Angle threshold
% % threshold_ptc_segment.search_method = "knn";
% % threshold_ptc_segment.search_value = 20;
% % 
% % 
% % threshold_voxel_segment.voxel_size = 0.05;
% % threshold_voxel_segment.angle = 15;
% % threshold_voxel_segment.residual = 0.025;
% % threshold_voxel_segment.distance = 0.025;
% % threshold_voxel_segment.min_ptc = 15;
% % threshold_voxel_segment.cosin_angle = cos(threshold_voxel_segment.angle*pi/180); %Angle threshold
% % threshold_voxel_segment.region_cosin_angle = cos(2.*threshold_voxel_segment.angle*pi/180); %Angle threshold
% % 
% % 
% % 
% % file_name = filename;
% 
% 
% % % Constant for pier
% % pier.component_min_area = 0.5;      % The minimum of an are of the component or surface = 1.0*0.5       
% % pier.section_min_edge_length = 0.5; % The minimum length of the edge of the structural component
% % pier.section_thickness = 2.0*threshold.sampling_step;      % Interval thcikness to extract the points for a cross-section identification
% % pier.section_distance = 0.25;          % A distance between two consecutive sections
% % pier.section_bbox_error = 0.2; % Deviation between cross-section representing by the bounding box
% % pier.segment_distance = 2.0*pier.section_distance; % distance between two consecutive sections on the plane
% % pier.column_min_length = 5.0*pier.section_distance; % The shortest column
% % pier.min_no_section = ceil(0.5*pier.column_min_length/pier.section_distance) ; % The minimum number of the cross-sections allowing to have deviation more than pier_threhold.section_error
% % 