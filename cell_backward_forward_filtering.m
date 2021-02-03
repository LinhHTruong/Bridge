function [Cells_Region] = cell_backward_forward_filtering(Tree, cells_plane, region_info, threshold, debug)
%% This function is to remove points on the boundary of the region out of the region
% because of over-segmentation
% 
% Input:

%       
% Output:
%       Region - cell ids[2x1], ptc_ids, surface_features(6x1)

% Developed by: Dr. Linh Truong-Hong, ORSG, Dept. GRS, TUDelft
%               
% %% Demo
% TREE = OQTR;
% Cells = Plane_Cells;
% THRESHOLD = THRESHOLD;
% region_info = Region_Info;


%% Estbalish the region structure
if debug
    tic
end
no_region = max(region_info(:,end));
Cells_Region = struct('cell',[], 'status', []);
% preallocate structure
for region_id = 1:no_region
    Cells_Region(region_id).status = [];
    Cells_Region(region_id).cell.id = [];
    
end
%% region_count = 1;
for i = 1:no_region
    % Retrieve region info
    mask = region_info(:,end) == i;
    region_cell_ids = region_info(mask,[1,2]);
    % Find cells on boundary of the region, which has less than 8 neibour
    % region
    for j = 1:size(region_cell_ids, 1)
          
        % Extract the points within the current cell
        region_cur_cell_id = region_cell_ids(j,:);
        
        % Searching neighbour cells
        neighbour_cell_ids = Query_Neighbour_Cells(Tree, region_cur_cell_id(1),region_info(:,1));
        % Extract the cells in/out region 
        mask = ismember(region_info(:,1),neighbour_cell_ids);
        neighbour_cell_peak_ids = region_info(mask,:);
        
        mask = ismember(neighbour_cell_peak_ids(:,end), i);
        in_region_cell_info = neighbour_cell_peak_ids(mask,:);
        out_region_cell_info = neighbour_cell_peak_ids(~mask,:);

        if isempty(out_region_cell_info) %Isolate region or on the boundary and no attached to any region 
            % Retrieve points in the peak shape of the current cell
            mask = ismember(cells_plane.cell_ids(:,[1,2]), region_cur_cell_id, 'rows');
            region_cur_cell_peak_ptc_ids = cells_plane.peak_info(mask).ptc_ids;
            
            if size(in_region_cell_info,1) == 8 %Interior cell -> no filtering
                % Cell peak features
                region_cur_cell_peak_features = cells_plane.peak_info(mask).peaks_features(1:end-1);
            else
                % Points in the current 
                region_cur_cell_peak_ptc_xyz = Tree.pts(region_cur_cell_peak_ptc_ids,1:3);

                % Retrieve the points in neighbour cells in the same region
                % -> fit a local surface -> use to filter the points within the current cell
                mask = ismember(cells_plane.cell_ids(:,[1,2]), in_region_cell_info(:,[1,2]), 'rows');
                in_region_neighbour_cell_peak_ptc_ids = vertcat(cells_plane.peak_info(mask).ptc_ids);
                in_region_neighbour_cell_peak_ptc_xyz = Tree.pts(in_region_neighbour_cell_peak_ptc_ids,1:3);
                
                [in_region_neighbour_cell_peak_cent, in_region_neighbour_cell_peak_normal] = wPCA('ptc', in_region_neighbour_cell_peak_ptc_xyz, 'angle_threshold', 0.5*pi/180, 'ini_normal_vector', []);
%                 in_region_neighbour_cell_peak_cent = mean(in_region_neighbour_cell_peak_ptc_xyz,1);
%                 in_region_neighbour_cell_peak_normal = eigenspace(in_region_neighbour_cell_peak_ptc_xyz,1);
                
                % Compute distance
                dist_region_ptc = abs(dist_3Dpoints_3Dplane(region_cur_cell_peak_ptc_xyz, [in_region_neighbour_cell_peak_cent, in_region_neighbour_cell_peak_normal]));
                
                % Filter inlier
                mask = dist_region_ptc <= threshold.max_distance;
                region_cur_cell_peak_ptc_ids = region_cur_cell_peak_ptc_ids(mask);
                if numel(region_cur_cell_peak_ptc_ids) >= threshold.min_num_pts
                    region_cur_cell_peak_ptc_xyz = region_cur_cell_peak_ptc_xyz(mask,:);
                    region_cur_cell_peak_features = [mean(region_cur_cell_peak_ptc_xyz,1), eigenspace(region_cur_cell_peak_ptc_xyz,1)];
                else
                    region_cur_cell_peak_features = inf;
                end
                clear in_region_neighbour_cell_peak_ptc_ids in_region_neighbour_cell_peak_ptc_xyz in_region_neighbour_cell_peak_center in_region_neighbour_cell_peak_normal
            end
            
            % Assign the new regions
            if isempty(Cells_Region(i).cell)
                count = 1;
            else
                count = size(vertcat(Cells_Region(i).cell.id),1) + 1;
            end
            
            % Update the region
            if ~isnan(region_cur_cell_peak_features) & ~isinf(region_cur_cell_peak_features)
                Cells_Region(i).status = 1;
                Cells_Region(i).cell(count).id = region_cur_cell_id;
                Cells_Region(i).cell(count).ptc_ids = region_cur_cell_peak_ptc_ids;
                Cells_Region(i).cell(count).surface_features = region_cur_cell_peak_features;
            end

            
        else % Either boundary (attached to another region) 
            % Retrieve points in the peak shape of the current cell
            mask = ismember(cells_plane.cell_ids(:,[1,2]), region_cur_cell_id, 'rows');
            region_cur_cell_peak_ptc_ids = cells_plane.peak_info(mask).ptc_ids;
            region_cur_cell_peak_ptc_xyz = Tree.pts(region_cur_cell_peak_ptc_ids,1:3);

            % Determine a local surface of the current region at the current cell, which is 
            % based on points within the cells'peaks in the same region: cells connected to the current cell
            if isempty(in_region_cell_info)
                mask = ismember(cells_plane.cell_ids(:, [1,2]), region_cur_cell_id, 'rows');
                region_local_surface_features = cells_plane.peak_info(mask).peaks_features(1:end-1);
            else
                % Searching interior cells
                region_neighbour_cell_peak_ptc_ids = region_cur_cell_peak_ptc_ids;
                for k = 1:numel(in_region_cell_info(:,1))
                    neighbour_cell_ids = Query_Neighbour_Cells(Tree, in_region_cell_info(k,1),region_cell_ids(:,1));
                    if numel(neighbour_cell_ids) ~= 8
                        % Filter the cells in the current region
                        mask = ismember(region_cell_ids(:,1),neighbour_cell_ids);
                        neighbour_cell_peak_ids = region_cell_ids(mask,:);
                        
                        % Retrieve the points in peak shape of neighbour
                        % cell
                        mask = ismember(cells_plane.cell_ids(:,[1,2]),neighbour_cell_peak_ids, 'rows');
                        neighbour_cell_peak_ptc_ids = vertcat(cells_plane.peak_info(mask).ptc_ids);
                        region_neighbour_cell_peak_ptc_ids = union(region_neighbour_cell_peak_ptc_ids,neighbour_cell_peak_ptc_ids); 
                        clear mask neighbour_cell_peak_ids neighbour_cell_peak_ptc_ids
                    end
                end
                % Compute features of the local surface of the current
                % region at the current cell
                region_neighbour_cell_peak_ptc_xyz = Tree.pts(region_neighbour_cell_peak_ptc_ids,1:3);
                local_surface_normal = eigenspace(region_neighbour_cell_peak_ptc_xyz, 1);
                region_local_surface_features = [mean(region_neighbour_cell_peak_ptc_xyz,1), local_surface_normal];
                clear local_surface_eigen_vectors local_surface_residual
            end
            % Compute distances of the points of the current cell to adjacent regions
            % Preallocate 
            out_region_ids = unique(out_region_cell_info(:,end));
            dist_region_ptc = zeros(numel(region_cur_cell_peak_ptc_ids),numel(out_region_ids) + 1);

            % For the local surface of the current region
            dist_region_ptc(:,1) = abs(dist_3Dpoints_3Dplane(region_cur_cell_peak_ptc_xyz, [region_local_surface_features(1:3), region_local_surface_features(4:6)]));

            % For the local surfaces of the neighbour region
            for k = 1:numel(out_region_ids)
                cur_out_region_id = out_region_ids(k);
                mask = out_region_cell_info(:,end) == cur_out_region_id;
                cur_out_region_cell_peak_ids = out_region_cell_info(mask,[1,2]);
                
                % Find all cells connected to 
                mask = region_info(:,end) == cur_out_region_id;
                cur_out_region_cell_info = region_info(mask,:);
                out_region_neighbour_cell_ids = Query_Neighbour_27Cells(Tree, cur_out_region_cell_peak_ids(:,1),cur_out_region_cell_info(:,1));
                % Retrieve cell ids + peak ids
                mask = ismember(cur_out_region_cell_peak_ids(:,1), out_region_neighbour_cell_ids);
                out_region_neighbour_cell_peak_ids = cur_out_region_cell_peak_ids(mask,[1,2]);
                
                % Compute features of the local surface of out_region and
                % distances from the points to the surface
                mask = ismember(cells_plane.cell_ids(:,[1,2]),out_region_neighbour_cell_peak_ids, 'rows');
                out_region_neighbour_cell_peak_ptc_ids = vertcat(cells_plane.peak_info(mask).ptc_ids);
                out_region_neighbour_cell_peak_ptc_xyz = Tree.pts(out_region_neighbour_cell_peak_ptc_ids,1:3);
                
                 [local_surface_cent, local_surface_normal] = wPCA('ptc', out_region_neighbour_cell_peak_ptc_xyz, 'angle_threshold', 0.5*pi/180, 'ini_normal_vector', []);
%              
                
%                 local_surface_normal = eigenspace(out_region_neighbour_cell_peak_ptc_xyz, 1);
                
                dist_region_ptc(:,k+1) = abs(dist_3Dpoints_3Dplane(region_cur_cell_peak_ptc_xyz, [local_surface_cent, local_surface_normal]));
                clear mask out_region_neighbour_cell_peak_ptc_ids out_region_neighbour_cell_peak_ptc_xyz local_surface_normal
            end
            
            % Filter the points for the regions
            all_region_ids = union(i, out_region_ids,'rows', 'stable');

            % Condition 1: no more threshold
            mask_1 = any((dist_region_ptc <= threshold.max_distance),2);
            % Condition 2: min distance
            [~, mask_2] = min(dist_region_ptc,[],2);
            mask = mask_1&mask_2;
            
            add_region_ptc_info = [region_cur_cell_peak_ptc_ids(mask), region_cur_cell_peak_ptc_xyz(mask,1:3), all_region_ids(mask_2(mask_1))];
            add_region_ids = unique(add_region_ptc_info(:,end));

            % Add cell and points to a new region when the number of the
            % points in the cell is larger than a predefined threshold         
            for k = 1:numel(add_region_ids)
                % Retrieve new points in a region
                add_region_id = add_region_ids(k);
                mask = add_region_ptc_info(:,end) == add_region_id;
                add_region_cell_ptc_ids = add_region_ptc_info(mask,1);
                if numel(add_region_cell_ptc_ids) >= threshold.min_num_pts
                    add_region_cell_ptc_xyz = add_region_ptc_info(mask,2:4);
                    if isempty(Cells_Region(add_region_id).cell)
                        count = 1;
                    else
                        count = size(vertcat(Cells_Region(add_region_id).cell.id),1) + 1;
                    end
                    % Update region
                    add_region_cell_normal = eigenspace(add_region_cell_ptc_xyz, 1);
                    Cells_Region(add_region_id).status = 1;
                    Cells_Region(add_region_id).cell(count).id = region_cur_cell_id;
                    Cells_Region(add_region_id).cell(count).ptc_ids = add_region_cell_ptc_ids;
                    Cells_Region(add_region_id).cell(count).surface_features = [mean(add_region_cell_ptc_xyz,1), add_region_cell_normal];
                end
                clear mask add_region_id add_region_cell_ptc_ids
            end
        end
    end
    %
end
if debug
    fprintf('Running time of the backward and forward filtering: %.2f seconds \n', toc);
end

%
% clear mask region_cell_ids cell_bound_info 
% clear interior_cell_ids count_j current_cell_id current_cell_order current_cell_ptc_ids current_cell_ptc_xyz
% clear add_neighbour_interior_cell_ids all_neighbour_interior_cell_ids bound_cell_ids
% clear i j k dist 
% clear in_region_current_cell_ptc_ids out_region_current_cell_ptc_ids in_region_current_cell_ptc_xyz out_region_current_cell_ptc_xyz
% clear no_region neighbour_cell_ids neighbour_interior_cell_ids neighbour_interior_cell_ptc_xyz region_count surface_center surface_normls      
