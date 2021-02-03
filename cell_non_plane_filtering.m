function [cells_plane, cells_region] = cell_non_plane_filtering(Tree, cells_plane, cells_region, threshold, debug)
%% This function is to examine if the points on non-plane peaks can be parts of the segment or not, 
% Which is based on the distance from the points to the local plane of the
% segments
% 
% Input:

%       
% Output:
%       Region - region id, cell ids, pts_ids, ptc_xyz
%       Cell - cell_ids, non_planar_ids, non_planar_xyz: not belong to
%       any segments
% Developed by: Dr. Linh Truong-Hong, ORSG, Dept. GRS, TUDelft
%               
% %% Demo
% TREE = OQTR;
% Cells = nonPlane_Cells;
% Cell_Region = Region;



   

%% Assign the points for defined regions
if debug
    tic
end
num_region = length(cells_region);
for region_id = 1:num_region
    
    % Extract the current region
    region_cell_peak_ids = vertcat(cells_region(region_id).cell.id);

    for j = 1:size(region_cell_peak_ids,1)
        % Searching neighbour cells out of the current region
        neighbour_cell_ids = Query_Neighbour_Cells(Tree, region_cell_peak_ids(j,1),cells_plane.cell_ids);

        if ~isempty(neighbour_cell_ids)
            
            % Retrieve surface features of the current cell in a current region
            region_cell_peak_pts = cells_region(region_id).cell(j).ptc_ids;
            region_cell_peak_xyz = Tree.pts(region_cell_peak_pts,1:3);
            region_cell_surface_features = cells_region(region_id).cell(j).surface_features;
            
            for k = 1:numel(neighbour_cell_ids)
                % Retrieve points in non-plane-cell
                neighbour_cell_id = neighbour_cell_ids(k);
                mask = ismember(cells_plane.cell_ids, neighbour_cell_id);
                non_plane_cell_ptc_ids = cells_plane.ptc_ids(mask).id;

                % Only consider if the number of cell larger than a predefined
                % threshold
                if (numel(non_plane_cell_ptc_ids) >= threshold.min_num_pts)
                    % Compute a distance from the points to the cell
                    non_plane_cell_ptc_xyz = Tree.pts(non_plane_cell_ptc_ids,1:3);
                    dist_ptc_2_cell = abs(dist_3Dpoints_3Dplane(non_plane_cell_ptc_xyz, region_cell_surface_features));

                    % Add new cell and points to the region: Condition 1: the
                    % number of the points
                     mask = dist_ptc_2_cell <= threshold.max_distance;
                     add_cell_ptc_ids = non_plane_cell_ptc_ids(mask);

                     if numel(add_cell_ptc_ids) >= threshold.min_num_pts
                        % Extract points to be added
                        add_cell_ptc_xyz = non_plane_cell_ptc_xyz(mask,:);

                        % Compute the distance between two cloud
                        min_dist_ptc_cell_2_cell = pdist2(region_cell_peak_xyz,add_cell_ptc_xyz);

                        % Condition 2: distance between two clouds
                        if min(min_dist_ptc_cell_2_cell(:)) <= 0.25*threshold.cell_size
                            % Compute surface feaures of the add cell
                            add_cell_ptc_normal = eigenspace(add_cell_ptc_xyz, 1);
                            if abs(cosine_vectors(region_cell_surface_features(4:6), add_cell_ptc_normal)) <= threshold.max_angle
                                add_cell_peak_ptc_xyz = union(region_cell_peak_xyz,add_cell_ptc_xyz, 'rows'); 
                                add_cell_ptc_normal = eigenspace(add_cell_peak_ptc_xyz, 1);
                            end
                            
                            % Update the add cell to the current region 
                            region_length = length(cells_region(region_id).cell);

                            cells_region(region_id).cell(region_length + 1).id = [neighbour_cell_id, 0];
                            cells_region(region_id).cell(region_length + 1).ptc_ids = add_cell_ptc_ids;
                            cells_region(region_id).cell(region_length + 1).surface_features = [mean(add_cell_ptc_xyz, 1), add_cell_ptc_normal];

                            % Update non-plane-cell: Cells
                            mask = ismember(cells_plane.cell_ids, neighbour_cell_id);
                            cells_plane.ptc_ids(mask).id = setdiff(non_plane_cell_ptc_ids, add_cell_ptc_ids);
                        end
                     end
                end
            end
        end
    end
end

 
% %% Update region: remove any empty region
% temp_Region = Cell_Region;
% Cell_Region = struct('id',[],'cell',[]);
% count = 1;
% for i=1:length(temp_Region)
%     if ~isempty(vertcat(temp_Region(i).cell.id))
%         Cell_Region(count).id = count;
%         Cell_Region(count).cell = temp_Region(i).cell;
%         count = count + 1;
%     end    
% end

if debug
    fprintf('Running time of the non_plane cell filtering: %.2f seconds \n', toc);
end
    %
% clear mask region_cell_ids cell_bound_info 
% clear interior_cell_ids count_j current_cell_id current_cell_order current_cell_ptc_ids current_cell_ptc_xyz
% clear add_neighbour_interior_cell_ids all_neighbour_interior_cell_ids bound_cell_ids
% clear i j k dist 
% clear in_region_current_cell_ptc_ids out_region_current_cell_ptc_ids in_region_current_cell_ptc_xyz out_region_current_cell_ptc_xyz
% clear no_region neighbour_cell_ids neighbour_interior_cell_ids neighbour_interior_cell_ptc_xyz region_count surface_center surface_normls      
