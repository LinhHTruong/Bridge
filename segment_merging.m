function new_cells_region = segment_merging(Tree, cells_region, threshold, debug)

%% This function is to merge the small segments together. At this step is only the planar segment
% Input:
% Output:
% Developed by: Dr. Linh Truong-Hong, ORSG, Dept. GRS, TUDelft
%               
% %% Demo
% Tree = OQTR;
% region = Region;
% threshold = THRESHOLD;
% debug = true;

% %% Summary the region
% num_region = length(region);
% region_stat = zeros(num_region,2);
% region_stat(:, 1) = [1:num_region]';
% for i = 1:num_region
%     region_stat(i,2) = size(vertcat(region(i).cell.id), 1);
% end
% [~, ids] = sort(region_stat(:, 2), "descend");
% Region_Ids = region_stat(ids, 1);
% clear num_region region_stat i ids
%% Compute features of each segment
region_features = cal_region_features(Tree, cells_region);

% Determine the region as a planar or not
mask = all(region_features.residual(:,2) <= 2.0*threshold.max_residual, 2);
plane_region_ids = region_features.ids(mask);
% curve_region_ids = Region_Features.ids(~mask);
% Re-order the region in its area
plane_region_area = region_features.bounding_box(plane_region_ids,1);
[~, sort_ids] = sort(plane_region_area, 'descend');
plane_region_ids = plane_region_ids(sort_ids);

clear mask sort_ids

%% Merging regions from different patches
if debug
    tic
end
% Preallocation
check_region_ids = plane_region_ids;
update_region = zeros(size(plane_region_ids,1),2);
update_region(:,1) = plane_region_ids;

% 
count = 1;
while ~isempty(check_region_ids)
%     numel(check_region_ids)
    if size(check_region_ids,1) == 1
        mask = ismember(update_region(:,1), check_region_ids, 'rows');
        update_region(mask, 2) = count;
        break
    end
    
    % Select initial region
    source_region_id = check_region_ids(1);
    % Update region 
    mask = ismember(update_region(:,1), source_region_id, 'rows');
    update_region(mask, 2) = count;
    
    % Update seeding region
    seed_region_ids = source_region_id;
    
    % Update check_patch_ind
    check_region_ids(1) = [];
    
    % Searching merging region
    while ~isempty(seed_region_ids)
        
        % Extract features of the ref_region
        cur_seed_region_id = seed_region_ids(1);
        
        % Determine a type of the seed region
        mask = ismember(plane_region_ids, cur_seed_region_id);

        if any(mask) %planar surface
            % Retrieve feature of the seed region
            mask = ismember(region_features.ids,cur_seed_region_id);
            cur_seed_region_cent = region_features.center(mask,:);
            cur_seed_region_normal = region_features.normal(mask,:);
            
            % Retrieve features of the check region
            target_region_ids = sort(check_region_ids);
            mask = ismember(region_features.ids, target_region_ids);
            target_region_cent = region_features.center(mask,:);
            target_region_normal = region_features.normal(mask,:);

            % Condition 1: Check similarity: direction -> normal deviation
            cosin_seed_check_regions = abs(cosine_vectors(cur_seed_region_normal, target_region_normal));
            mask = cosin_seed_check_regions >= cos(deg2rad(2.0*threshold.max_angle));
            target_region_cent = target_region_cent(mask,:);
            target_region_ids = target_region_ids(mask,:);

            % Check orthogonal distance between regions
            if ~isempty(target_region_ids)
                
                % Condition 2: Compute distance
                dist_seed_check_regions = dist_3Dpoints_3Dplane(target_region_cent, [cur_seed_region_cent, cur_seed_region_normal]);
                mask = abs(dist_seed_check_regions) <= 0.5*threshold.min_dist_2_planes;
                target_region_ids = target_region_ids(mask);
                % Condition 3: Check continuity
                if numel(target_region_ids) > 0
                    % Determine boundary cells of the current_seed
                    cur_seed_region_cell_ids = vertcat(cells_region(cur_seed_region_id).cell.id);
                    
                    % check the current region connect to the target
                    for j = 1:numel(target_region_ids)
                        % Extract the bounds in current_target_region
                        cur_target_region_id = target_region_ids(j);
                        cur_target_region_cell_ids = vertcat(cells_region(cur_target_region_id).cell.id);
                        
                        % Find neighbour cells of bound cells of the
                        % current seed region located in the target
                        seed_target_neighbour_cell_ids = Query_Neighbour_Cells(Tree, cur_seed_region_cell_ids(:,1), cur_target_region_cell_ids(:,1));
                        target_seed_neighbour_cell_ids = Query_Neighbour_Cells(Tree, cur_target_region_cell_ids(:,1), cur_seed_region_cell_ids(:,1));

                        % Overlap between adjacent cells from
                        % a seeding region vs. ones from a target region
                        if (~isempty(seed_target_neighbour_cell_ids))&&(~isempty(target_seed_neighbour_cell_ids))
                            % Extract points in boundaries area
                            % From current seed region
                            mask = ismember(cur_seed_region_cell_ids, target_seed_neighbour_cell_ids);
                            cur_seed_region_bound_cell_ptc_ids = vertcat(cells_region(cur_seed_region_id).cell(mask).ptc_ids);
                            cur_seed_region_bound_cell_ptc_xyz = Tree.pts(cur_seed_region_bound_cell_ptc_ids,1:3);
                            
                            % From current target region
                            mask = ismember(cur_target_region_cell_ids, seed_target_neighbour_cell_ids);
                            cur_target_region_bound_cell_ptc_ids = vertcat(cells_region(cur_target_region_id).cell(mask).ptc_ids);
                            cur_target_region_bound_cell_ptc_xyz = Tree.pts(cur_target_region_bound_cell_ptc_ids,1:3);
                            
                            % Compute distance of join points to its
                            % fitting plane
                            bound_ptc = union(cur_seed_region_bound_cell_ptc_xyz,cur_target_region_bound_cell_ptc_xyz, "rows");
                            [wptc_center, wptc_normal] = wPCA('ptc', bound_ptc);
                            dist_bound_ptc_plane = dist_3Dpoints_3Dplane(bound_ptc, [wptc_center, wptc_normal]);
                            residual_bound_ptc_plane = sqrt(sum(dist_bound_ptc_plane.^2)/numel(dist_bound_ptc_plane));
                            
                            % Determine threshold
                            mask = ismember(region_features.ids, [cur_seed_region_id,cur_target_region_id]);
                            threshold_residual = region_features.residual(mask,:);
                            if (mean(abs(dist_bound_ptc_plane)) <= max(threshold_residual(:,1)))&&(residual_bound_ptc_plane <= max(threshold_residual(:,2)))
                                % Merging - Update new region id
                                mask = ismember(update_region(:,1), cur_target_region_id, 'rows');
                                update_region(mask, 2) = count;
                                % Update seeding region
                                seed_region_ids = union(seed_region_ids, cur_target_region_id);
                                % Remove target region to be merge
                                check_region_ids = setdiff(check_region_ids,cur_target_region_id, 'stable');
                            end
                        end
                    end
                end
            end
            
        end
        % remove the seed region to be check
        seed_region_ids = setdiff(seed_region_ids, cur_seed_region_id);
    end  
    count = count + 1;
end
fprintf('Running time for merging the small segments: %.2f seconds \n', toc);
clear all_region_cell_ids
%% Update 
tic
new_cells_region = struct('cell',[], 'status', []);
for i = 1:length(unique(update_region(:,2)))
    mask = update_region(:,2) == i;
    
    % Retrive information from current re
    old_region_ids = update_region(mask,1);
    % Update a new region
    new_cells_region(i).status = cells_region(old_region_ids(1)).status;
    new_cells_region(i).cell = cells_region(old_region_ids(1)).cell;
    if numel(old_region_ids) > 1  
        for j = 2:numel(old_region_ids)
            % Retrieve cell ids in the new region
            new_region_cell_ids = vertcat(new_cells_region(i).cell.id);
            
            % Retrieve cell ids in old regions to check any duplication
            old_region_cell_ids = vertcat(cells_region(old_region_ids(j)).cell.id);
            for k = 1:size(old_region_cell_ids,1)
                % Check if any duplication
                mask = ismember(new_region_cell_ids, old_region_cell_ids(k,:), 'rows');
                if any(mask)
                    % Retrieve points and features in the cell in the new region
                    new_region_cell_ptc_ids = new_cells_region(i).cell(mask).ptc_ids;
                    
                    % Retrieve points and features in the cell in the old region
                    old_region_cell_ptc_ids = cells_region(old_region_ids(j)).cell(k).ptc_ids;
                    
                    % Update new data for the cell
                    cell_update_ptc_ids = union(new_region_cell_ptc_ids,old_region_cell_ptc_ids);
                    cell_update_ptc_xyz = Tree.pts(cell_update_ptc_ids,1:3);
                    cell_center = mean(cell_update_ptc_xyz);
                    cell_normal = eigenspace(cell_update_ptc_xyz, 1);

                    % Update the points in the old region to the new region
                    new_cells_region(i).cell(mask).ptc_ids = cell_update_ptc_ids;
                   
                    new_cells_region(i).cell(mask).surface_features = [cell_center, cell_normal];

                else
                    
                    % Update the new cell to the region
                    count_id = length(new_cells_region(i).cell);
                    new_cells_region(i).cell(count_id + 1).id = old_region_cell_ids(k,:);
                    new_cells_region(i).cell(count_id + 1).ptc_ids = cells_region(old_region_ids(j)).cell(k).ptc_ids;
                    new_cells_region(i).cell(count_id + 1).surface_features = cells_region(old_region_ids(j)).cell(k).surface_features;
                end
            end
        end
    end    
end
if debug
    fprintf('Running time for updating a new segment: %.2f seconds \n', toc);
end

