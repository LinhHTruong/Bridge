function [global_region_info, cells_plane] = cell_region_growing_segmentation(Tree, cells_plane, threshold, option, debug)
% This function is to segmentation a point cloud which particularly focus on
% the bridge deck
%   - Region growing based quadtree or cells, in which deviation ofthe 
% normal vector and heigh difference are considered as the conditions for the
% growing progress while the residual is the threshold for selecting
% seeding cells

% Input:
% Output:
% Developed by: Dr. Linh Truong-Hong, ORSG, Dept. GRS, TUDelft
%               
% %% Demo
% Tree = OQTR;
% cells = Plane_Cells;
% threshold = THRESHOLD;
% option = 'top';

%% Extract the cells
if debug
    tic
end
% %%
mask = cells_plane.cell_ids(:,end) == 1;
cells_ids = unique(cells_plane.cell_ids(mask,1));
cells_peaks_features = [];
for i = 1:numel(cells_ids)
    cell_id = cells_ids(i);
    mask = (cells_plane.cell_ids(:,1) == cell_id)&(cells_plane.cell_ids(:,end) == 1);
    cell_peaks_ids = cells_plane.cell_ids(mask,:);
    peaks_features = vertcat(cells_plane.peak_info(mask).peaks_features);
    
    % Remove any cell_peaks nearly vertical
    cosin_cells_peaks_oz = abs(cosine_vectors(peaks_features(:,4:6), threshold.nz));
    mask = cosin_cells_peaks_oz >= cos(deg2rad(threshold.hor_plane_max_angle));
    cell_peaks_ids = cell_peaks_ids(mask,:);
         
    % Filtering the peaks can consider as seeding cell peak
    if ~isempty(cell_peaks_ids)
        % Update cell_peak features
        peaks_features = peaks_features(mask,:);
        
        if size(cell_peaks_ids,1) == 1
            cells_peaks_features = [cells_peaks_features;[cell_peaks_ids(:,[1,2]),peaks_features]]; 
        else
            if strcmp(option, 'bottom')
                % Start from the bottom peak
                [~, id] = max(cell_peaks_ids(:,2));
            else 
                % Start from the top peak
                [~, id] = min(cell_peaks_ids(:,2));
            end  
            cells_peaks_features = [cells_peaks_features;[cell_peaks_ids(id,[1,2]),peaks_features(id,:)]]; 
        end
    end

end
clear cells_ids cell_id mask i id cell_peaks_ids peaks_features

%% Initial region
% Extract initial seeding cells
ini_seed_cells = cells_peaks_features(:,[1,end]);

% Global region
global_region_info = inf(size(cells_peaks_features,1),3);
global_region_info(:,1:2) = cells_peaks_features(:,[1,2]);

checked_cells = zeros(size(cells_peaks_features,1),2);
checked_cells(:,1) = cells_peaks_features(:,1);
checked_cells(:,2) = cells_peaks_features(:,end);

region_no = 0;
%
while ~isempty(ini_seed_cells)
    % Initial seeding voxel
    [~,mask] = min(ini_seed_cells(:,2));
    seed_cell_id = ini_seed_cells(mask,1);
    % Seeding region
    seed_region_ids = seed_cell_id;
    seed_region_ids = seed_region_ids(:);
    % Current region
    cur_region = seed_cell_id;
    cur_region = cur_region(:);
    % Update checked voxel
    mask = ismember(checked_cells(:,1),seed_cell_id); 
    checked_cells(mask,:) = [];
    % Region growing
    while ~isempty(seed_region_ids) 
        
        cur_seed_cell_id = seed_region_ids(1);
        % Searching neighbour
        neighbour_cell_ids = Query_Neighbour_Cells(Tree,cur_seed_cell_id,checked_cells(:,1));
        
        % Check if neighbour_cell_ids in checked_cells
        mask = ismember(neighbour_cell_ids, checked_cells);
        neighbour_cell_ids = neighbour_cell_ids(mask);
        if ~isempty(neighbour_cell_ids)
            % Extract current and neighbour cells information
            mask = ismember(cells_peaks_features(:,1), cur_seed_cell_id);
            cur_seed_cell_cent = cells_peaks_features(mask,3:5); 
            cur_seed_cell_normal = cells_peaks_features(mask,6:8);

            % Retrieve features of neighbour cells
            mask = ismember(cells_peaks_features(:,1), neighbour_cell_ids);
            neighbour_cell_cent = cells_peaks_features(mask,3:5); 
            neighbour_cell_normal = cells_peaks_features(mask,6:8);
           
            % Check conditions
            cosin_seed_neighbour_cells = abs(cosine_vectors(cur_seed_cell_normal, neighbour_cell_normal));
            dist_seed_neighbour_cells = abs(dist_3Dpoints_3Dplane(neighbour_cell_cent, [cur_seed_cell_cent, cur_seed_cell_normal]));
            mask = (abs(cosin_seed_neighbour_cells) > cos(deg2rad(threshold.max_angle)))&...
                   (dist_seed_neighbour_cells <= threshold.max_distance);
            add_cell_ids = neighbour_cell_ids(mask);
            
            % Add cells to current region
            if numel(add_cell_ids) > 0
                % Update a current region
                cur_region = union(cur_region, add_cell_ids); 
                % Remove new add_seeding out of the checked voxel
                mask = ismember(checked_cells(:,1), add_cell_ids);
                checked_cells(mask,:) = [];
                
                % Check if the cell can be add to seeding region
                mask = ismember(cells_peaks_features(:,1), add_cell_ids);
                add_deck_cell_res = cells_peaks_features(mask,end);
          
                % Update a seeding region
                mask = add_deck_cell_res <= threshold.max_residual;
                if any(mask)
                    add_seed_cell_ids = add_cell_ids(mask);
                    seed_region_ids = union(seed_region_ids,add_seed_cell_ids);
                end
            end
        end
        % Remove current_seeding_voxel_ind
        seed_region_ids = setdiff(seed_region_ids, cur_seed_cell_id);
    end 
    
    % Update initial seeding 
    mask = ismember(ini_seed_cells(:,1), cur_region);
    ini_seed_cells(mask,:) = [];
    
    % Update the current region
    if threshold.min_num_cell <= numel(cur_region)
        mask = ismember(global_region_info(:,1), cur_region);
        global_region_info(mask,3) = region_no + 1;
        region_no = region_no + 1;
    end
end

clear num_cell region_no checked_cell seed_cell_id seed_region cur_region
clear cur_seed_cell_id neighbour_cell_id cur_seed_cell_normal 
clear neighbour_cell_normal add_cell_id  
clear cosin_seed_neighbour_cells add_seed_cell_ind mask 
%% Update region ids
% Remove non-segment cells
mask = isinf(global_region_info(:,end));
global_region_info(mask,:) = [];

num_region = max(global_region_info(:,end));
region_count = histcounts(global_region_info(:,end), num_region);
region_stat(:, 1) = unique(global_region_info(:,end));
region_stat(:, 2) = region_count;

[~, ids] = sort(region_stat(:, 2), "descend");
region_stat = region_stat(ids,:);

% Update region_info
temp_region_ids = inf(size(global_region_info, 1),1);
for i = 1:size(region_stat,1)
    mask = ismember(global_region_info(:, end), region_stat(i,1));
    temp_region_ids(mask) = i;
end
%
global_region_info(:,end) = temp_region_ids;
clear num_region region_count region_stat ids mask temp_region_ids

if debug
    fprintf('Running time of the segmentation: %.2f seconds \n', toc);
end


