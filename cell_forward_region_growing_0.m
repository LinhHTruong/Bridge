function [cells, cell_region] = cell_forward_region_growing(Tree, cells, cell_region, threshold, debug)
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
% %% Demo
% Tree = OQTR;
% Cells = Plane_Cells;
% Cell_Region = Region;
% threshold = THRESHOLD;
%% Extract the region cell peak ids
if strcmpi(debug, 'true')
    tic
end
no_region = length(cell_region);
regions_cells_peaks_info = [];
for region_id = 1:no_region
    cell_peak_ids = vertcat(cell_region(region_id).cell.id);    
    cell_peak_ids(:,3) = region_id;
    % Cummulative region
    regions_cells_peaks_info = [regions_cells_peaks_info;cell_peak_ids];
end 
% Update the staus of the peak
mask = ismember(cells.cell_ids(:,[1,2]), regions_cells_peaks_info(:,[1,2]), 'rows');
cells.cell_ids(mask, end) = 0;
non_regions_cells_peaks_info = cells.cell_ids(~mask, :);
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
tic
num_region = length(cell_region);
for region_id = 1:num_region
    % Assign non-segment cells
    check_non_regions_cells_peaks_info = non_regions_cells_peaks_info;
    
    % Extract the region
    region_cell_peak_ids = vertcat(cell_region(region_id).cell.id);
    check_cell_peak_ids = region_cell_peak_ids;
%    
    flag = true;
    while flag

%         initial_check_cell_peak_ids = check_cell_peak_ids;
        for count = 1:size(check_cell_peak_ids,1)
            % Current cell of the current region
            check_cell_peak_id = check_cell_peak_ids(1,:);

            % Searching neighbour cells out of the current region
            neighbour_cell_ids = Query_Edge_Neighbour_Cells(Tree, check_cell_peak_id(1),check_non_regions_cells_peaks_info(:,1));

            % Check each peaks of the cells
            if ~isempty(neighbour_cell_ids)
                % Retrieve cell_peak
                mask = ismember(check_non_regions_cells_peaks_info(:,1),neighbour_cell_ids);
                neighbour_cell_peak_ids = check_non_regions_cells_peaks_info(mask,:);
                mask = neighbour_cell_peak_ids(:,end) == 1;
                neighbour_cell_peak_ids = neighbour_cell_peak_ids(mask,:);
%                 flag_neighbour_cell_peak_ids(mask) = true;
               
                % Retrieve surface features of the current cell
                mask = ismember(region_cell_peak_ids, check_cell_peak_id, 'rows');
                check_cell_peak_ptc_ids = cell_region(region_id).cell(mask).ptc_ids;
                check_cell_peak_ptc_xyz = Tree.pts(check_cell_peak_ptc_ids,1:3);
                check_cell_surf_features = cell_region(region_id).cell(mask).surface_features;

                % Check each peak of the cell
%                 flag_add_neighbour_cell_peak_ids = false(size(neighbour_cell_peak_ids,1),1);
                for i = 1:size(neighbour_cell_peak_ids,1)
                    % Retrieve the points within the peak
                    cur_neighbour_cell_peak_id = neighbour_cell_peak_ids(i,[1,2]);
                    mask = ismember(cells.cell_ids(:,[1,2]), cur_neighbour_cell_peak_id(:,[1,2]), 'rows');
                    cur_neighbour_cell_peak_ptc_ids = cells.peak_info(mask).ptc_ids;
                    cur_neighbour_cell_peak_ptc_xyz = Tree.pts(cur_neighbour_cell_peak_ptc_ids,1:3);
                    cur_neighbour_cell_surf_features = cells.peak_info(mask).peaks_features;

                    % Calculate a distance from the points to the current cell
                    dist_neighbour_cell_ptc_cell = abs(dist_3Dpoints_3Dplane(cur_neighbour_cell_peak_ptc_xyz, check_cell_surf_features));

                    % Add cell and points
                    mask = dist_neighbour_cell_ptc_cell <= threshold.max_distance;
                    add_cell_peak_ptc_ids = cur_neighbour_cell_peak_ptc_ids(mask);

                    if (numel(add_cell_peak_ptc_ids) >= threshold.min_num_pts)

                        add_cell_peak_ptc_xyz = cur_neighbour_cell_peak_ptc_xyz(mask,:);

                        % Compute a distance between points in a check cell and ones in add-cell
                        threshold_range_dist = 0.5*norm(check_cell_surf_features(1:3) - cur_neighbour_cell_surf_features(1:3));%0.25*THRESHOLD.cell_size)
                        min_dist_check_add_cell_ptc = pdist2(add_cell_peak_ptc_xyz, check_cell_peak_ptc_xyz);
                        mask = any(min_dist_check_add_cell_ptc <= threshold_range_dist, 2);
                        
                        
                        if (sum(mask) >= threshold.min_num_pts)
                            
                            % Mark the neighbour cell peak id was added
%                             flag_add_neighbour_cell_peak_ids(i) = true;
                            
                            % Check if the cell in the region already
                            mask = ismember(vertcat(cell_region(region_id).cell.id), cur_neighbour_cell_peak_id,'rows');
                            if find(mask == 1)
                                cur_region_cell_ptc_ids = cell_region(region_id).cell(mask).ptc_ids;
                                cur_region_cell_ptc_xyz = Tree.pts(cur_region_cell_ptc_ids,1:3);

                                % Update new data for the cell
                                cur_region_cell_update_ptc_ids = union(cur_region_cell_ptc_ids,add_cell_peak_ptc_ids);
                                cur_region_cell_update_ptc_xyz = union(cur_region_cell_ptc_xyz,add_cell_peak_ptc_xyz, 'rows');

                                % Update cell region
                                cell_region(region_id).cell(mask).ptc_ids = cur_region_cell_update_ptc_ids;
                                cell_region(region_id).cell(mask).surface_features = [mean(cur_region_cell_update_ptc_xyz, 1), eigenspace(cur_region_cell_update_ptc_xyz, 1)];

                            else
                                 
                                % Check if the adding cell containing the small number of the points, which lead to
                                % sigularity (the normal is large different from the check cell: IF no, use both data points
                                add_cell_peak_normal = eigenspace(add_cell_peak_ptc_xyz, 1);
                                if abs(cosine_vectors(check_cell_surf_features(4:6), add_cell_peak_normal)) <= cos(deg2rad(threshold.max_angle))
                                    add_cell_peak_ptc_xyz = union(check_cell_peak_ptc_xyz, add_cell_peak_ptc_xyz, 'rows'); 
                                    add_cell_peak_normal = eigenspace(add_cell_peak_ptc_xyz, 1);
                                end

                                % Update the cell region 
                                region_length = length(cell_region(region_id).cell);
                                cell_region(region_id).cell(region_length + 1).id = cur_neighbour_cell_peak_id;
                                cell_region(region_id).cell(region_length + 1).ptc_ids = add_cell_peak_ptc_ids;
                                cell_region(region_id).cell(region_length + 1).surface_features = [mean(add_cell_peak_ptc_xyz, 1), add_cell_peak_normal];

                                % Update a new adding cell to check_cell_peak_id and region_cell_peak_ids
                                cosin_check_neighour_cells = abs(cosine_vectors(check_cell_surf_features(4:6),cur_neighbour_cell_surf_features(4:6)));
                                if all(~isinf(cur_neighbour_cell_surf_features)) & (cosin_check_neighour_cells >= cos(deg2rad(threshold.hor_plane_max_angle)))
                                    check_cell_peak_ids = union(check_cell_peak_ids, cur_neighbour_cell_peak_id,'rows','stable');
                                end
                                region_cell_peak_ids = union(region_cell_peak_ids, cur_neighbour_cell_peak_id,'rows','stable');
                            end

                            % Update ptc_ids in the Cells because part of the points were added to a new region
                            mask = ismember(cells.cell_ids(:,[1,2]), cur_neighbour_cell_peak_id, 'rows');
                            % Update cell_peak status: 0 - inactive
                            cells.cell_ids(mask,end) = 0;

                            % Update points in the peak shape
                            non_region_cell_peak_ptc_ids = cells.peak_info(mask).ptc_ids;
                            non_region_cell_peak_remain_ptc_ids = setdiff(non_region_cell_peak_ptc_ids,add_cell_peak_ptc_ids);
                            cells.peak_info(mask).ptc_ids = non_region_cell_peak_remain_ptc_ids;

                        end
                        
%                     else
                        % no more peek to be check for a case top -
                        % down
                    end
%                     end
                end
                
%                 % Update cell_peak ids were added for 
%                 add_neighbour_cell_peak_ids = neighbour_cell_peak_ids(flag_add_neighbour_cell_peak_ids,:);
%                 mask = ismember(check_non_regions_cells_peaks_info, add_neighbour_cell_peak_ids, 'rows');
%                 flag_neighbour_cell_peak_ids(mask) = true;

            end
            
%             % Remove check_cell_peak_id
            check_cell_peak_ids = setdiff(check_cell_peak_ids, check_cell_peak_id, 'rows','stable');
            
        end
        
        if isempty(check_cell_peak_ids) || isempty(check_non_regions_cells_peaks_info)
            flag = false;
        end
            
    end
end
if strcmpi(debug, 'true')
    fprintf('Running time of the forward growing segmentation: %.2f seconds \n', toc);
end
% clear num_region region_cell_id region_cell_prop index 
% clear bound_cell_id bound_cell_surface_centroid bound_cell_surface_normls 
% clear neighbour_bound_cell_id neighbour_bound_cell_ptc_id neighbour_bound_cell_ptc_xyz
% clear dist cell_add_ptc_xyz mask nodeLeaf cell_all_ptc_ids cell_all_ptc_xyz wptc_center wptc_normal
