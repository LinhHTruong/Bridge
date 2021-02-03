function region_info = plane_outlier_removal(Tree, region_info, region_ids, threshold)
%% This function is to remove oulier points of the region.
% The points of the region are considered as outlier points if they are away from the intersection line
% 
% Input:
%       Tree                    : Data structure stored entire data
%       cells_regions           : data structure stored cells region
%       regions_ids             : a pair of the region intersection
%
% Output:
%
% Demo:
% region_info = cells_region;
% region_info = ptc_segment_info;
% region_ids = comp_link_list(2,[1,2]);
%% Retrieve the points of the region
% Extract the points indicies
if isstruct(region_info)
    region_1_ptc_ids = vertcat(region_info(region_ids(1)).cell.ptc_ids);
    region_2_ptc_ids = vertcat(region_info(region_ids(2)).cell.ptc_ids);
else
    mask = region_info(:,2) == region_ids(1);
    region_1_ptc_ids = region_info(mask,1);
    mask = region_info(:,2) == region_ids(2);
    region_2_ptc_ids = region_info(mask,1);   
end

% Retrieve points
region_1_ptc_xyz = Tree.pts(region_1_ptc_ids,1:3);
region_2_ptc_xyz = Tree.pts(region_2_ptc_ids,1:3);

% Determine an intersection between two surface
region_1_cent = mean(region_1_ptc_xyz, 1);
region_1_normal = eigenspace(region_1_ptc_xyz,1);
region_2_cent = mean(region_2_ptc_xyz, 1);
region_2_normal = eigenspace(region_2_ptc_xyz,1);
[intersect_line, ~] = plane_plane_intersection([region_1_cent, region_1_normal], [region_2_cent, region_2_normal], threshold.max_angle);

% Gathering the normal
region_surface = [region_1_cent, region_1_normal;...
                  region_2_cent, region_2_normal];
clear region_1_ptc_ids  region_1_ptc_xyz region_2_ptc_ids region_2_ptc_xyz 
%% Extract the outlier points
for i = 1:numel(region_ids)
    % Retrieve the points of the region
%     i = 1;
    if isstruct(region_info)
        region_cell_ids = vertcat(region_info(region_ids(i)).cell.id);
        [region_bound_cell_ids, ~] = get_boundary_cells(Tree, region_cell_ids(:,1));
        mask = ismember(region_cell_ids(:,1), region_bound_cell_ids);
        region_bound_ptc_ids = vertcat(region_info(region_ids(i)).cell(mask).ptc_ids);
    else
        mask = region_info(:,2) == region_ids(i);
        region_bound_ptc_ids = region_info(mask,1); %Use all points of the segment
    end
    region_bound_ptc_xyz = Tree.pts(region_bound_ptc_ids,1:3);

    % Determine the line 
    % Project points to the fitting surface
    region_bound_proj_ptc_xyz = proj_3Dpoints_3Dplane(region_bound_ptc_xyz, region_surface(i,:));
    
    % Project the center to the line
    region_proj_center = proj_point_on_3Dline(region_surface(i,1:3), intersect_line);
    
    % Cal sign distances
    d_ptc_line = dist_3Dpoints_3Dline(region_bound_proj_ptc_xyz, intersect_line, region_proj_center - region_surface(i,1:3));
    
    % Classify the points based on sign distance
    sign_dist = sign(d_ptc_line);
    sign_cluster_ids = unique(sign_dist);
    cluster_count = histcounts(sign_dist, numel(sign_cluster_ids));

    % Determine inlier and outlier ptc
    [~, max_id] = max(cluster_count);
    sign_inlier_cluster_ids = sign_cluster_ids(max_id);
    mask = sign_dist == sign_inlier_cluster_ids;
%     inlier_ptc_ids = ptc_ids(mask);
    region_outier_ptc_ids = region_bound_ptc_ids(~mask);
    
    % Update inlier points
    if isstruct(region_info)
        % Update the cells
        for j = 1:size(region_bound_cell_ids,1)
            % Retrieve the points within the bound cell
            region_cell_id = region_bound_cell_ids(j,:);
            mask = ismember(region_cell_ids(:,1), region_cell_id);
            local_id = find(mask);
            region_ptc_ids = region_info(region_ids(i)).cell(local_id).ptc_ids;

            % Remove outlier points
            mask = ismember(region_ptc_ids, region_outier_ptc_ids);
            if any(mask)
                region_info(region_ids(i)).cell(local_id).ptc_ids = region_ptc_ids(~mask);
            end
        end
    else
        % update array
        mask = ismember(region_info(:,1), region_outier_ptc_ids);
        region_info(mask,2) = 0;
    end
end
    
% 
% %% Write txt file
% 
% outputFile = 'segment 1.txt';
% write_txt(outputFile, region_1_bound_ptc_xyz)
% 
% outputFile = 'segment 2.txt';
% write_txt(outputFile, region_2_bound_ptc_xyz)
% 
% outputFile = 'outlier segment 1.txt';
% write_txt(outputFile, region_bound_ptc_xyz(~mask,:))
% 
% region_1_ptc_ids = vertcat(cells_region(region_ids(1)).cell.ptc_ids);
% region_1_bound_ptc_xyz = Tree.pts(region_1_ptc_ids,1:3);
% 
% outputFile = 'after segment 1.txt';
% write_txt(outputFile, region_1_bound_ptc_xyz)
%     


