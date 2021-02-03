function [cells_plane, cells_region] = cell_boundary_filtering(Tree, cells_plane, cells_region, threshold, debug)
%% This function is to filter points on the boundary of the region which may be incorrect asignment to the region
% 
% Input:

%       
% Output:
%       Cells - cell ids[2x1], ptc_ids, surface_features(6x1)
%       Cell_Region

% Developed by: Dr. Linh Truong-Hong, ORSG, Dept. GRS, TUDelft
%               
% %% Demo
% TREE = OQTR;
% Cells = Plane_Cells;
% Cell_Region = Region;
% THRESHOLD = THRESHOLD;


%% Establish a list cell peak region
if debug
    tic
end
regions_cells_info = [];
for region_id = 1:length(cells_region)
    cell_peak_ids = vertcat(cells_region(region_id).cell.id);    
    cell_peak_ids(:,3) = region_id;
    % Cummulative region
    regions_cells_info = [regions_cells_info;cell_peak_ids];
end 
%% Extract all boundary cells
no_region = length(cells_region);
boundary_cells = [];
for region_id = 1:no_region
    % Retrieve region info
    mask = regions_cells_info(:,end) == region_id;
    region_connecting_bound_cell_peak_ids = regions_cells_info(mask,:);
    % Find cells on boundary of the region, which has less than 8 neibour
    % region
    for j = 1:size(region_connecting_bound_cell_peak_ids, 1)
        region_current_cell_id = region_connecting_bound_cell_peak_ids(j,:);
        
        % Searching neighbour cells
        neighbour_cell_ids = Query_Neighbour_Cells(Tree, region_current_cell_id(1),region_connecting_bound_cell_peak_ids(:,1));
        if numel(neighbour_cell_ids) < 8
            boundary_cells = [boundary_cells;region_current_cell_id];
        end
    end
end
%% Extract sharing cell_peaks on the boundary of the region: sharing by at least two regions
[~,ia, ib] = unique(regions_cells_info(:,[1,2]),'rows','stable');
[count_cell_peak, ~, count_ids] = histcounts(ib,numel(ia));
bound_cell_peak_ids = count_cell_peak(count_ids)>1;
regions_sharing_bound_cells_info = regions_cells_info(bound_cell_peak_ids, :);

%% Extract connecting cell_peaks on the boundary of the region: sharing by at least two regions
mask = ismember(boundary_cells, regions_sharing_bound_cells_info, 'rows');
regions_connecting_bound_cells_info = boundary_cells(~mask,:);

%% Filtering for connecting cells
no_region = length(cells_region);
for region_id = 1:no_region
    % Retrieve cell_peak id
    region_cell_peak_ids = vertcat(cells_region(region_id).cell.id);
    % Retrieve region info
    mask = regions_connecting_bound_cells_info(:,end) == region_id;
    region_connecting_bound_cell_peak_ids = regions_connecting_bound_cells_info(mask,[1,2]);
    % Find cells on boundary of the region, which has less than 8 neibour
    % region
    for j = 1:size(region_connecting_bound_cell_peak_ids, 1)
          
        % Extract the points within the current cell
        region_current_cell_id = region_connecting_bound_cell_peak_ids(j,:);
        
        % Searching neighbour cells
        neighbour_cell_ids = Query_Neighbour_Cells(Tree, region_current_cell_id(1),regions_cells_info(:,1));
        
        % Extract the cells in/out region 
        mask = ismember(regions_cells_info(:,1),neighbour_cell_ids);
        neighbour_cell_peak_ids = regions_cells_info(mask,:);
        mask = ismember(neighbour_cell_peak_ids(:,end), region_id);
        in_region_cell_info = neighbour_cell_peak_ids(mask,:);
        out_region_cell_info = neighbour_cell_peak_ids(~mask,:);

        if ~isempty(out_region_cell_info)  
            % Retrieve points in of the current cell in the current region
            mask = ismember(region_cell_peak_ids(:,[1,2]), region_current_cell_id, 'rows');
            region_current_bound_cell_peak_ptc_ids = cells_region(region_id).cell(mask).ptc_ids;
            region_current_bound_cell_peak_ptc_xyz = Tree.pts(region_current_bound_cell_peak_ptc_ids,1:3);

            % Determine a local surface of the current region at the current cell, which is 
            % based on points within the cells'peaks in the same region: cells connected to the current cell
            if isempty(in_region_cell_info)
                region_local_surface_features = cells_region(region_id).cell(mask).surface_features;
            else
                % Compute features of the local surface of the current
                % region at the current cell based on its neighbour cells
                mask = ismember(region_cell_peak_ids(:,[1,2]), in_region_cell_info(:,[1,2]), 'rows');
                in_region_neighbour_cell_peak_ptc_ids = cells_region(region_id).cell(mask).ptc_ids;
                in_region_neighbour_cell_peak_ptc_xyz = Tree.pts(in_region_neighbour_cell_peak_ptc_ids,1:3);

                local_surface_normal = eigenspace(in_region_neighbour_cell_peak_ptc_xyz, 1);
                region_local_surface_features = [mean(in_region_neighbour_cell_peak_ptc_xyz,1), local_surface_normal];
                clear mask local_surface_normal in_region_neighbour_cell_peak_ptc_ids in_region_neighbour_cell_peak_ptc_xyz
            end

            % Compute distances of the points of the current cell to adjacent regions
            % Preallocate 
            out_region_ids = unique(out_region_cell_info(:,end));
            dist_region_ptc = zeros(numel(region_current_bound_cell_peak_ptc_ids),numel(out_region_ids) + 1);

            % For the local surface of the current region
            dist_region_ptc(:,1) = abs(dist_3Dpoints_3Dplane(region_current_bound_cell_peak_ptc_xyz, [region_local_surface_features(1:3), region_local_surface_features(4:6)]));

            % For the local surfaces of the neighbour region
            for k = 1:numel(out_region_ids)
                % Retrieve cell peaks in the same neighbour region
                current_out_region_id = out_region_ids(k);
                mask = out_region_cell_info(:,end) == current_out_region_id;
                current_out_region_cell_peak_ids = out_region_cell_info(mask,[1,2]);
                
                % Find all cells connected to 
                out_region_cell_peak_ids = vertcat(cells_region(current_out_region_id).cell.id);
                out_region_neighbour_cell_ids = Query_Neighbour_27Cells(Tree, current_out_region_cell_peak_ids(:,1),out_region_cell_peak_ids(:,1));
                % Retrieve cell ids + peak ids
                mask = ismember(out_region_cell_peak_ids(:,1), out_region_neighbour_cell_ids);
                out_region_neighbour_cell_peak_ids = out_region_cell_peak_ids(mask,[1,2]);
                
                % Compute features of the local surface of out_region and
                % distances from the points to the surface
                mask = ismember(out_region_cell_peak_ids(:,[1,2]),out_region_neighbour_cell_peak_ids, 'rows');
                out_region_neighbour_cell_peak_ptc_ids = vertcat(cells_region(current_out_region_id).cell(mask).ptc_ids);
                out_region_neighbour_cell_peak_ptc_xyz = Tree.pts(out_region_neighbour_cell_peak_ptc_ids,1:3);
                local_surface_normal = eigenspace(out_region_neighbour_cell_peak_ptc_xyz, 1);
                
                % Check angle between the neighbour cell peaks to the current
                cosin_in_out_regions = abs(cosine_vectors(region_local_surface_features(4:6), local_surface_normal));
                if cosin_in_out_regions >= threshold.max_angle
                    dist_region_ptc(:,k+1) = abs(dist_3Dpoints_3Dplane(region_current_bound_cell_peak_ptc_xyz, [mean(out_region_neighbour_cell_peak_ptc_xyz,1), local_surface_normal]));
                else
                    dist_region_ptc(:,k+1) = 1000; % No longer consider this cell_peaks
                end
                clear mask out_region_neighbour_cell_peak_ptc_ids out_region_neighbour_cell_peak_ptc_xyz local_surface_normal
            end
            
            % Filter the points for the regions
            all_region_ids = union(region_id, out_region_ids,'rows', 'stable');

            % Condition 1: no more threshold
            mask_1 = any((dist_region_ptc <= threshold.max_distance),2);
            % Condition 2: min distance
            [~, mask_2] = min(dist_region_ptc,[],2);
            mask = mask_1&mask_2;
            
            % Determine points to be added to the region
            add_region_cell_peak_ptc_info = [region_current_bound_cell_peak_ptc_ids(mask), all_region_ids(mask_2(mask_1))];
            add_region_ids = unique(add_region_cell_peak_ptc_info(:,end));
            region_cell_peak_fragment_ptc_ids = region_current_bound_cell_peak_ptc_ids(~mask);
            
            % Add cell and points to a new region when the number of the
            % points in the cell is larger than a predefined threshold         
            for k = 1:numel(add_region_ids)
                % Retrieve new points in a region
                add_region_id = add_region_ids(k);
                mask = add_region_cell_peak_ptc_info(:,end) == add_region_id;
                add_region_cell_peak_ptc_ids = add_region_cell_peak_ptc_info(mask,1);
                if numel(add_region_cell_peak_ptc_ids) >= threshold.min_num_pts
                    % Retrieve adding points 
                    add_region_cell_peak_ptc_xyz = Tree.pts(add_region_cell_peak_ptc_ids,1:3);
                    
                    % Update data for cell
                    if add_region_id == region_id
                        % Update for a current region
                        if ~isequal(region_current_bound_cell_peak_ptc_ids,add_region_cell_peak_ptc_ids) 
                            % Need to update
                            mask = ismember(region_cell_peak_ids(:,[1,2]),region_current_cell_id, 'rows');
                            cells_region(region_id).cell(mask).ptc_ids = add_region_cell_peak_ptc_ids;
                            cells_region(region_id).cell(mask).surface_features = [mean(add_region_cell_peak_ptc_xyz,1), eigenspace(add_region_cell_peak_ptc_xyz,1)];
                        end
                    else
                        % Update for neighbour region
                        out_region_cell_peak_ids = vertcat(cells_region(add_region_id).cell.id);
                        
                        % Check if the cell-peak already in the region 
                        mask = ismember(out_region_cell_peak_ids(:,[1,2]),region_current_cell_id, 'rows');
                        if any(mask)
                            % The cell peak is already in the region - Update points
                            out_region_current_cell_peak_ptc_ids = cells_region(region_id).cell(mask).ptc_ids;
                            out_region_current_cell_peak_update_ptc_ids = union(out_region_current_cell_peak_ptc_ids, add_region_cell_peak_ptc_ids);
                            out_region_current_cell_peak_update_ptc_xyz = Tree.pts(out_region_current_cell_peak_update_ptc_ids,1:3);

                            % Update cell_peak
                            cells_region(add_region_id).cell(mask).ptc_ids = out_region_current_cell_peak_update_ptc_ids;
                            cells_region(add_region_id).cell(mask).surface_features = [mean(out_region_current_cell_peak_update_ptc_xyz,1), eigenspace(out_region_current_cell_peak_update_ptc_xyz,1)];
                            clear out_region_current_cell_peak_ptc_ids out_region_current_cell_peak_update_ptc_ids out_region_current_cell_peak_update_ptc_xyz
                        else
                            % The cell peak is NOT already in the region
                            region_length = length(cells_region(add_region_id).cell);

                            cells_region(add_region_id).cell(region_length + 1).id = region_current_cell_id;
                            cells_region(add_region_id).cell(region_length + 1).ptc_ids = add_region_cell_peak_ptc_ids;
                            cells_region(add_region_id).cell(region_length + 1).surface_features = [mean(add_region_cell_peak_ptc_xyz, 1), eigenspace(add_region_cell_peak_ptc_xyz,1)];
                            clear region_length add_region_cell_ptc_xyz
                        end
                    end

                else
                    % Points no add to the region although they own distances smaller than a threshold
                    region_cell_peak_fragment_ptc_ids = [region_cell_peak_fragment_ptc_ids;add_region_cell_peak_ptc_ids];
                end
                clear mask add_region_cell_ptc_ids
            end
            
            % Return the points within the current cell-peak are not
            % assigned to the region back to the original cell
            mask = ismember(cells_plane.cell_ids(:,[1,2]), region_current_cell_id, 'rows');
            
            current_cell_ptc_ids = cells_plane.peak_info(mask).ptc_ids;
            cells_plane.peak_info(mask).ptc_ids = unique([current_cell_ptc_ids;region_cell_peak_fragment_ptc_ids]);
            clear mask current_cell_ptc_ids fragment_region_ids

        end
    end
    %
end

fprintf('Running time of the connecting cell filtering: %.2f seconds \n', toc);

%% Filtering for sharing cells
tic
for i = 1:size(regions_sharing_bound_cells_info,1)
    % Retrieve a current sharing_bound_cell
    cell_peak_info = regions_sharing_bound_cells_info(i,:);
    
    % Searching neighbour cells
    neighbour_cell_ids = Query_Neighbour_Cells(Tree, cell_peak_info(1),regions_cells_info(:,1));
    mask = ismember(regions_cells_info(:,1), neighbour_cell_ids);
    neighbour_cell_peak_info = regions_cells_info(mask,:);
    % Determine the regions occupied the cell
    sharing_regions_ids = unique(neighbour_cell_peak_info(:,end));
    
    % Extract points within the cell_peak assigned to the region
    sharing_regions_cell_peak_ptc_ids = [];
    flag_region = false(numel(sharing_regions_ids), 1);
    for j = 1:numel(sharing_regions_ids)
        sharing_region_id = sharing_regions_ids(j);
        mask = ismember(vertcat(cells_region(sharing_region_id).cell.id), cell_peak_info(1:2), 'rows');
        if any(mask)
            flag_region(j) = true;
            sharing_region_cell_peak_ptc_ids = cells_region(sharing_region_id).cell(mask).ptc_ids;
            sharing_regions_cell_peak_ptc_ids = [sharing_regions_cell_peak_ptc_ids;sharing_region_cell_peak_ptc_ids];
        end
    end
    clear sharing_region_id mask sharing_region_cell_peak_ptc_ids
    
    % Filtering the regions to be considered
    sharing_regions_ids = sharing_regions_ids(flag_region);
    mask = ismember(neighbour_cell_peak_info(:,end),sharing_regions_ids);
    neighbour_cell_peak_info = neighbour_cell_peak_info(mask,:);
    clear flag_region

    % Compute distance from the points to the region
    sharing_regions_cell_peak_ptc_xyz = Tree.pts(sharing_regions_cell_peak_ptc_ids,1:3);
    dist_region_ptc = zeros(numel(sharing_regions_cell_peak_ptc_ids),numel(sharing_regions_ids));

    for j = 1:numel(sharing_regions_ids)
        % Retrieve connect cell-peaks of the region to the current-sharing-bound-region
        sharing_region_id = sharing_regions_ids(j);
        mask = ismember(neighbour_cell_peak_info(:,end),sharing_region_id);
        sharing_cell_peak_ids = neighbour_cell_peak_info(mask,[1,2]);

        % Retrieve points within those cells
        mask = ismember(vertcat(cells_region(sharing_region_id).cell.id), sharing_cell_peak_ids(:,1:2), 'rows');
        sharing_cell_peak_ptc_ids = vertcat(cells_region(sharing_region_id).cell(mask).ptc_ids);
        sharing_cell_peak_ptc_xyz = Tree.pts(sharing_cell_peak_ptc_ids,1:3);
        
        % Compute a local surface
        [sharing_cell_peak_surface_center, sharing_cell_peak_surface_normal] = wPCA('ptc', sharing_cell_peak_ptc_xyz, 'ini_normal_vector', [0,0,1]);
         
        % Compute distances
        dist_region_ptc(:,j) = abs(dist_3Dpoints_3Dplane(sharing_regions_cell_peak_ptc_xyz, [sharing_cell_peak_surface_center, sharing_cell_peak_surface_normal]));
           
        clear sharing_region_id mask sharing_cell_peak_ptc_ids sharing_cell_peak_ptc_xyz sharing_cell_peak_surface_center sharing_cell_peak_surface_normal
        
    end

    % Re-determine regions for the points
    % Condition 1: no more threshold
    mask_1 = any((dist_region_ptc <= threshold.max_distance),2);
    % Condition 2: min distance
    [~, mask_2] = min(dist_region_ptc,[],2);
    mask = mask_1&mask_2;

    % Determine points to be added to the region
    adjust_regions_cell_peak_ptc_info = [sharing_regions_cell_peak_ptc_ids(mask), sharing_regions_ids(mask_2(mask_1))];
    adjust_region_ids = unique(adjust_regions_cell_peak_ptc_info(:,end));
    cell_peak_fragment_ptc_ids = sharing_regions_cell_peak_ptc_ids(~mask); % return to original

    
    % Add cell and points to a new region when the number of the
    % points in the cell is larger than a predefined threshold         
    for k = 1:numel(adjust_region_ids)
        % Retrieve new points in a region
        adjust_region_id = adjust_region_ids(k);
        mask = adjust_regions_cell_peak_ptc_info(:,end) == adjust_region_id;
        adjust_region_cell_peak_ptc_ids = adjust_regions_cell_peak_ptc_info(mask,1);
        if numel(adjust_region_cell_peak_ptc_ids) >= threshold.min_num_pts
            % Retrieve adding points 
            adjust_region_cell_peak_ptc_xyz = Tree.pts(adjust_region_cell_peak_ptc_ids,1:3);
            
            % Compute surface features
            adjust_region_cell_peak_surface_center = mean(adjust_region_cell_peak_ptc_xyz,1);
            adjust_region_cell_peak_surface_normal = eigenspace(adjust_region_cell_peak_ptc_xyz,1);    
            
            % Update data for cell
            mask = ismember(vertcat(cells_region(adjust_region_id).cell.id), cell_peak_info(:,1:2), 'rows');
            cells_region(adjust_region_id).cell(mask).ptc_ids = adjust_region_cell_peak_ptc_ids;
            cells_region(adjust_region_id).cell(mask).surface_features = [adjust_region_cell_peak_surface_center, adjust_region_cell_peak_surface_normal];
            
            clear adjust_region_cell_peak_ptc_xyz adjust_region_cell_peak_surface_center adjust_region_cell_peak_surface_normal  
        
        else
            % Points no add to the region although they own distances smaller than a threshold
            cell_peak_fragment_ptc_ids = [cell_peak_fragment_ptc_ids;adjust_region_cell_peak_ptc_ids];
        end
        clear mask adjust_region_id  adjust_region_cell_peak_ptc_ids
    end

    % Return the points within the current cell-peak are not
    % assigned to the region back to the original cell
    mask = ismember(cells_plane.cell_ids(:,[1,2]), cell_peak_info(:,1:2), 'rows');
    current_cell_ptc_ids = cells_plane.peak_info(mask).ptc_ids;
    cells_plane.peak_info(mask).ptc_ids = unique([current_cell_ptc_ids;cell_peak_fragment_ptc_ids]);
    clear mask current_cell_ptc_ids cell_peak_fragment_ptc_ids
end
if debug
    fprintf('Running time of the sharing cell filtering: %.2f seconds \n', toc);
end
