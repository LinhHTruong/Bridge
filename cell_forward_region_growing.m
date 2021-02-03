function [cells_plane, cells_region] = cell_forward_region_growing(Tree, cells_plane, cells_region, threshold, debug)
%% This function is to re-obtain the points of the current region that locate in peak shape of the cells which is not segments
% Note: 
%   - Only cells are non-segment considered
%   - This step allow the points on the boundary can be located within multiple segments
% Input:
%       TREE        - Data structure of the cell
%       Cell        - data structure of cells with the plane: Only Cell.ptc_ids and Cell.ptc_xyz are used
%       Cell_Region - data structure of cells were assined in the region:  region id, cell ids, pts_ids, ptc_xyz
%       THRESHOLD   - Thresholds
%       
% Output:
%       Cell        - data structure of cells with the plane: if the cell
%       assigned to the region, the points in the surface were removed
%       Cell_Region - New cell id associated to the surface and points
%       within surface were added
%
% Developed by: Dr. Linh Truong-Hong, ORSG, Dept. GRS, TUDelft
%               
%% Demo
% Tree = OQTR;
% cells_plane = cells_plane;
% cell_region = Region;
% threshold = THRESHOLD;
% debug = false
%% Extract the region cell peak ids
if debug
    tic
end
no_region = length(cells_region);
regions_cells_peaks_info = [];
for region_id = 1:no_region
    cell_peak_ids = vertcat(cells_region(region_id).cell.id);    
    cell_peak_ids(:,3) = region_id;
    % Cummulative region
    regions_cells_peaks_info = [regions_cells_peaks_info;cell_peak_ids];
end 
% Update the staus of the peak
mask = ismember(cells_plane.cell_ids(:,[1,2]), regions_cells_peaks_info(:,[1,2]), 'rows');
cells_plane.cell_ids(mask, end) = 0;
non_regions_cells_peaks_info = cells_plane.cell_ids(~mask, :);
clear region_id cell_peak_ids mask
%% Extract all boundary cells
regions_bound_cells_peaks_info = [];
for region_id = 1:no_region
    % Retrieve region info
    mask = regions_cells_peaks_info(:,end) == region_id;
    region_cell_peak_ids = regions_cells_peaks_info(mask,:);
    % Find cells on boundary of the region, which has less than 8 neibour
    % region
    for j = 1:size(region_cell_peak_ids, 1)
        region_cell_peak_id = region_cell_peak_ids(j,:);
        % Searching neighbour cells
        neighbour_cell_ids = Query_Neighbour_Cells(Tree, region_cell_peak_id(1),region_cell_peak_ids(:,1));
        if numel(neighbour_cell_ids) < 8
            regions_bound_cells_peaks_info = [regions_bound_cells_peaks_info;region_cell_peak_id];
        end
    end
end
clear mask region_cell_peak_ids region_cell_peak_id j neighbour_cell_ids
%% Remove peaks of interior cells
mask = ismember(regions_cells_peaks_info, regions_bound_cells_peaks_info, 'rows');
regions_interior_cells_peaks_info = regions_cells_peaks_info(~mask,:);
mask = ismember(non_regions_cells_peaks_info(:,1), regions_interior_cells_peaks_info(:,1), 'rows');
non_regions_cells_peaks_info(mask,:) = [];
clear mask regions_interior_cells_peaks_info

%% Searching cells out off the region containing points fitting the local surface of the region and then growing

num_region = length(cells_region);
for region_id = 1:num_region
    % Assign non-segment cells

    % Extract the region
    check_cell_peak_ids = vertcat(cells_region(region_id).cell.id);
    check_non_regions_cells_peaks_info = non_regions_cells_peaks_info;
    
    flag = true;
    while flag
        
        % Filtering cell_peaks
        mask = ismember(check_non_regions_cells_peaks_info(:,1),check_cell_peak_ids(:,1)); 
        check_non_regions_cells_peaks_info = check_non_regions_cells_peaks_info(~mask,:);
        check_non_regions_cells_ids = unique(check_non_regions_cells_peaks_info(:,1));
        add_check_non_regions_cells = [];
        flag_check_non_regions_cells = false(size(check_non_regions_cells_peaks_info,1),1);
        %
        for count = 1:size(check_cell_peak_ids,1)
            % Current cell of the current region
            check_cell_peak_id = check_cell_peak_ids(count,:);
            check_cell_bounds = Tree.cell_bounds(check_cell_peak_id(1),:);
            check_cell_poly_xy = check_cell_bounds([1,2; 4,2; 4,5; 1,5; 1,2]);
            

            % Searching neighbour cells out of the current region
            neighbour_cell_ids = Query_Edge_Neighbour_Cells(Tree, check_cell_peak_id(1), check_non_regions_cells_ids);

            if isempty(neighbour_cell_ids)
%                 check_cell_peak_ids = setdiff(check_cell_peak_ids, check_cell_peak_id, 'rows','stable');
                continue;
            end
            % Retrieve cell_peak
            mask = ismember(check_non_regions_cells_peaks_info(:,1),neighbour_cell_ids);
            flag_check_non_regions_cells(mask) = true;
            neighbour_cell_peak_ids = check_non_regions_cells_peaks_info(mask,:);
            mask = neighbour_cell_peak_ids(:,end) == 1;
            neighbour_cell_peak_ids = neighbour_cell_peak_ids(mask,:);
            
            % Check the points within the neighbour cell can be merge to the current region
%             flag_neighbour_cell = false(size(neighbour_cell_peak_ids,1),1);
            for i = 1:size(neighbour_cell_peak_ids,1)
                % Find cells in the region connect to neighbour cells
                neighbour_check_cell_ids = Query_Edge_Neighbour_Cells(Tree, neighbour_cell_peak_ids(i,1), check_cell_peak_ids(:,1));
                mask = ismember(check_cell_peak_ids(:,1), neighbour_check_cell_ids);
                neighbour_check_cell_peak_ids = check_cell_peak_ids(mask,:);

                % Retrieve the points within the neighbour_cell_peak_ids
                mask = ismember(cells_plane.cell_ids(:,[1,2]), neighbour_cell_peak_ids(i,[1,2]), 'rows');
                neighbour_cell_peak_ptc_ids = cells_plane.peak_info(mask).ptc_ids;
                neighbour_cell_peak_ptc_xyz = Tree.pts(neighbour_cell_peak_ptc_ids,1:3);
                
                % Check if it is necessary to find the inlier point or not
                if numel(neighbour_cell_peak_ptc_ids) < threshold.min_num_pts
                    continue;
                end
                
                % Estimate a local surface from the neighbour_check_cell
                mask = ismember(vertcat(cells_region(region_id).cell.id), neighbour_check_cell_peak_ids, 'rows');
                neighbour_check_cell_peak_ptc_ids = vertcat(cells_region(region_id).cell(mask).ptc_ids);
                neighbour_check_cell_peak_ptc_xyz = Tree.pts(neighbour_check_cell_peak_ptc_ids,1:3);
                if size(neighbour_check_cell_peak_ids,1) == 1
                    neighbour_check_cell_peak_surf = cells_region(region_id).cell(mask).surface_features;
                else
                    neighbour_check_cell_peak_ptc_cent = mean(neighbour_check_cell_peak_ptc_xyz,1);
                    neighbour_check_cell_peak_ptc_normal = eigenspace(neighbour_check_cell_peak_ptc_xyz,1);
                    neighbour_check_cell_peak_surf = [neighbour_check_cell_peak_ptc_cent,neighbour_check_cell_peak_ptc_normal];
                    clear neighbour_check_cell_peak_ptc_cent neighbour_check_cell_peak_ptc_normal  
                end
                % Check if any point fit to the local surface                
                dist_neighbour_ptc_2_local_surf = dist_3Dpoints_3Dplane(neighbour_cell_peak_ptc_xyz(:,1:3), neighbour_check_cell_peak_surf);
                mask = abs(dist_neighbour_ptc_2_local_surf) <= threshold.max_distance;
                cell_peak_add_ptc_ids = neighbour_cell_peak_ptc_ids(mask);
                cell_peak_add_ptc_xyz = neighbour_cell_peak_ptc_xyz(mask,:);
                
                % Check the short distance between two clouds
                cell_peak_add_proj_ptc_xyz = proj_3Dpoints_3Dplane(cell_peak_add_ptc_xyz, neighbour_check_cell_peak_surf);
%                 size(neighbour_check_cell_peak_ptc_xyz)
%                 size(cell_peak_add_proj_ptc_xyz)
                dist_ptc_2_ptc = pdist2(neighbour_check_cell_peak_ptc_xyz,cell_peak_add_proj_ptc_xyz);
                if all(dist_ptc_2_ptc(:) > threshold.cell_size/3)
                    continue;
                end
                clear cell_peak_add_proj_ptc_xyz dist_ptc_2_ptc neighbour_check_cell_peak_ptc_ids neighbour_check_cell_peak_ptc_xyz
%                 % Iterative searching points
%                 flag_add_ptc = false(numel(neighbour_cell_peak_ptc_ids),1);
%                 flag_select_ptc = false(numel(neighbour_cell_peak_ptc_ids),1);
%                 flag = true;
% %                 ini_num_ptc = 0;
%                 offset_val = threshold.cell_size/4;
%                 
%                 while flag & (offset_val <= threshold.cell_size)
%                     % Searching the point within the bounds
%                     check_cell_poly_xy_offset = polygon_offset(check_cell_poly_xy, offset_val);
%                     [in_ids,on_ids] = inpolygon(neighbour_cell_peak_ptc_xyz(:,1),neighbour_cell_peak_ptc_xyz(:,2),check_cell_poly_xy_offset(:,1),check_cell_poly_xy_offset(:,2));
%                     mask = in_ids|on_ids;
%                     cell_poly_ptc_ids = neighbour_cell_peak_ptc_ids(mask);
%                     cell_poly_ptc_xyz = neighbour_cell_peak_ptc_xyz(mask,:);
%                     
%                     % Remove the points have been selected
%                     mask = ismember(cell_poly_ptc_ids, neighbour_cell_peak_ptc_ids(~flag_select_ptc));
%                     cell_poly_ptc_ids = cell_poly_ptc_ids(mask);
%                     cell_poly_ptc_xyz = cell_poly_ptc_xyz(mask,:);
% 
%                     if numel(cell_poly_ptc_ids)  < 0.5*threshold.min_num_pts
%                         flag = false;
%                     else
%                     
%                         % Check if the points are inlier
%                         dist_neighbour_ptc_2_local_surf = dist_3Dpoints_3Dplane(cell_poly_ptc_xyz(:,1:3), neighbour_check_cell_peak_surf);
%                         mask = abs(dist_neighbour_ptc_2_local_surf) <= threshold.max_distance;
%                         cell_poly_add_ptc_ids = cell_poly_ptc_ids(mask);
%                         
% %                         cell_poly_add_ptc_xyz = cell_poly_ptc_xyz(mask,:);
%                         clear dist_neighbour_ptc_2_local_surf
%                         % Set new polygon for next searching
%                         if numel(cell_poly_add_ptc_ids) >= 0.5*threshold.min_num_pts
%                             % Mark inlier points
%                             mask = ismember(neighbour_cell_peak_ptc_ids,cell_poly_add_ptc_ids);
%                             flag_add_ptc(mask) = true;
% 
%                             % Update points for next iteration
%                             offset_val = offset_val + offset_val;
%                             
%                         else
%                             flag = false;
%                         end
%                         clear cell_poly_add_ptc_ids cell_poly_add_ptc_xyz
%                     end
%                 end
                
                % Check if the points can be added
%                 cell_ptc_ids = neighbour_cell_peak_ptc_ids;
%                 cell_ptc_xyz = neighbour_cell_peak_ptc_xyz;
                
%                 % Iterative searching points
%                 flag_cell_peak_ptc = false(numel(neighbour_cell_peak_ptc_ids),1);
%                 flag = true;
%                 offset_val = max(0.075,threshold.cell_size/3);
%                 while flag
%                     % Searching the point within the bounds
%                     check_cell_poly_xy_offset = polygon_offset(check_cell_poly_xy, offset_val);
%                     [in_ids,on_ids] = inpolygon(cell_ptc_xyz(:,1),cell_ptc_xyz(:,2),check_cell_poly_xy_offset(:,1),check_cell_poly_xy_offset(:,2));
%                     mask = in_ids|on_ids;
%                     cell_poly_ptc_ids = cell_ptc_ids(mask);
%                     cell_poly_ptc_xyz = cell_ptc_xyz(mask,:);
%                     
%                     if numel(cell_poly_ptc_ids) < 0.5*threshold.min_num_pts
%                         flag = false;
%                     else
%                     
%                         % Update the points in the cell
%                         cell_ptc_ids = cell_ptc_ids(~mask);
%                         cell_ptc_xyz = cell_ptc_xyz(~mask,:);
%                         clear in_ids on_ids mask
% 
%                         % Check if the points are inlier
%                         dist_neighbour_ptc_2_local_surf = dist_3Dpoints_3Dplane(cell_poly_ptc_xyz(:,1:3), neighbour_check_cell_peak_surf);
%                         mask = abs(dist_neighbour_ptc_2_local_surf) <= threshold.max_distance;
%                         cell_poly_add_ptc_ids = cell_poly_ptc_ids(mask);
%                         cell_poly_add_ptc_xyz = cell_poly_ptc_xyz(mask,:);
%                         clear dist_neighbour_ptc_2_local_surf
%                         % Set new polygon for next searching
%                         if numel(cell_poly_add_ptc_ids) >= 0.5*threshold.min_num_pts
%                             % Mark inlier points
%                             mask = ismember(neighbour_cell_peak_ptc_ids,cell_poly_add_ptc_ids);
%                             flag_cell_peak_ptc(mask) = true;
% 
%                             % New polygon
%                             new_bounds = [min(cell_poly_add_ptc_xyz,[], 1),max(cell_poly_add_ptc_xyz,[], 1)];
%                             check_cell_poly_xy = new_bounds([1,2; 4,2; 4,5; 1,5; 1,2]);
%                             clear new_bounds
%                         else
%                             flag = false;
%                         end
%                         clear cell_poly_add_ptc_ids cell_poly_add_ptc_xyz
%                     end
%                 end
                
                % Retrieve add_cell
%                 cell_peak_add_ptc_ids = neighbour_cell_peak_ptc_ids(flag_add_ptc);
%                 cell_peak_add_ptc_xyz = neighbour_cell_peak_ptc_xyz(flag_add_ptc,:);
                
                if numel(cell_peak_add_ptc_ids) >= threshold.min_num_pts

%                       flag_neighbour_cell(i) = true;
                    % Prevent sigularity
                    cell_peak_add_ptc_normal = eigenspace(cell_peak_add_ptc_xyz, 1);
                    if abs(cosine_vectors(neighbour_check_cell_peak_surf(4:6), cell_peak_add_ptc_normal)) <= cos(deg2rad(threshold.max_angle)) 
                        cell_peak_add_ptc_normal = neighbour_check_cell_peak_surf(4:6);
                    else
                        % Update the seeding cell for next integration
                        add_check_non_regions_cells = [add_check_non_regions_cells;neighbour_cell_peak_ids(i,:)]; 
%                   
                    end
                    % Add to a current region 
                    region_length = length(cells_region(region_id).cell);
                    cells_region(region_id).cell(region_length + 1).id = neighbour_cell_peak_ids(i,[1,2]);
                    cells_region(region_id).cell(region_length + 1).ptc_ids = cell_peak_add_ptc_ids;
                    cells_region(region_id).cell(region_length + 1).surface_features = [mean(cell_peak_add_ptc_xyz, 1), cell_peak_add_ptc_normal];

                    % Update ptc_ids in the Cells because part of the points were added to a new region
                    mask = ismember(cells_plane.cell_ids(:,[1,2]), neighbour_cell_peak_ids(i,[1,2]), 'rows');
                    cells_plane.cell_ids(mask,end) = 0; % 0 - deactive
                end

            end

            % Update check_non_regions_cells_peaks_info
%             if any(flag_neighbour_cell)
                remove_non_region_cell_peak = neighbour_cell_peak_ids;%(flag_neighbour_cell,:);
                mask = ismember(check_non_regions_cells_peaks_info(:,[1,2]),remove_non_region_cell_peak(:,[1,2]), 'rows');
                check_non_regions_cells_peaks_info(mask,:) = [];
                check_non_regions_cells_ids = unique(check_non_regions_cells_peaks_info(:,1));
%             end
        end
        
        %% Remove check_cell_peak_id  
        if isempty(add_check_non_regions_cells) || isempty(check_non_regions_cells_ids)
            flag = false;
        else
            check_cell_peak_ids = add_check_non_regions_cells(:,[1,2]);
        end
    end
end
if debug
    fprintf('Running time of the forward growing segmentation: %.2f seconds \n', toc);
end
% clear num_region region_cell_id region_cell_prop index 
% clear bound_cell_id bound_cell_surface_centroid bound_cell_surface_normls 
% clear neighbour_bound_cell_id neighbour_bound_cell_ptc_id neighbour_bound_cell_ptc_xyz
% clear dist cell_add_ptc_xyz mask nodeLeaf cell_all_ptc_ids cell_all_ptc_xyz wptc_center wptc_normal
