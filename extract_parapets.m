function supstruct = extract_parapets(Tree, supstruct, road_surface, struct_threshold, threshold)
%% This function is to extract the parapets (vehicle and pedestrian)
% Input:
%     Tree                - Data structure stored 2D cells
%     super_structure     - Data structure of the bridge stored a final result
%     struct_threshold    - Threshold related to the structure
%     threshold           - Threshold related to the data processing

% Output:
%     super_structure
% Demo
% Tree = OQTR;
% supstruct = SuperStructure;
% road_surface = Road_Surface;
% threshold = THRESHOLD;
%% Constant
min_num_cells = struct_threshold.bridge_min_length/threshold.cell_size;
bin_width_dist = struct_threshold.footpath_min_width/3;
%% Retrieve leaf nodes
leaf_cell_ids = Node_Leaf(Tree);
leaf_cell_bounds = Tree.cell_bounds(leaf_cell_ids,:);
leaf_cell_cent = 0.5*(leaf_cell_bounds(:,[4,5,6]) + leaf_cell_bounds(:,[1,2,3]));

%% Extract indices of the road and foot path

% Road surface
% mask = cellfun(@(s) contains('Road Surface', s), supstruct(1).Description);
% road_surf_id = find(mask);
% Foot path
mask = cellfun(@(s) contains(s, 'Footpath'), supstruct(1).Description);
footpath_ids = find(mask);

%% Retrieve the road surface
% road_surf_pts_ids = vertcat(supstruct(1).Component(road_surf_id).cell.ptc_ids);
% road_surf_pts_xyz = Tree.pts(road_surf_pts_ids,1:3);
% % road_surface_cent = mean(road_surface_ptc_xyz,1);
% % road_surface_normal = eigenspace(road_surface_ptc_xyz,1);
% % road_surface_tangent = eigenspace(road_surface_ptc_xyz,2);
% road_mbb = min2DBoundingBox(road_surf_pts_xyz(:,1:2)');
% road_surface_cent = [mean(road_mbb.vertices,1), 0];
% road_surface_normal = cross([road_mbb.long_edge_vector, 0], [road_mbb.short_edge_vector,0]);
% road_surface_tangent = [road_mbb.long_edge_vector,0];

clear road_surface_pts_ids road_surface_pts_xyz road_mbb

%% Extract the road curb

if ~isempty(footpath_ids)
    for i = 1: numel(footpath_ids)
        
        %% Retrieve cells of the footpath
        footpath_id = footpath_ids(i);
        footpath_cell_ids = vertcat(supstruct(1).Component(footpath_id).cell.ids);
        footpath_cell_pts_ids = vertcat(supstruct(1).Component(footpath_id).cell.ptc_ids);
        footpath_cell_pts_xyz = Tree.pts(footpath_cell_pts_ids,1:3);
        
        % Extract cells within the footpath
        footpath_mmb = min2DBoundingBox(footpath_cell_pts_xyz(:,1:2)');
        footpath_poly = [footpath_mmb.vertices; footpath_mmb.vertices(1,:)];
        [in_ids, on_ids] = inpolygon(leaf_cell_cent(:,1), leaf_cell_cent(:,2), footpath_poly(:,1), footpath_poly(:,2));
        parapet_cell_ids = leaf_cell_ids(in_ids | on_ids); 
        clear footpath_cell_pts_ids footpath_cell_pts_xyz footpath_mmb footpath_poly in_ids on_ids
        
        % Determine the points within the cell and located above the footpath
        parapet = struct('cell',[]);
        no_parapet_cell = 1;
        for j = 1:numel(parapet_cell_ids)
            % Retrieve the parapet_cell
            parapet_cell_id = parapet_cell_ids(j);
            parapet_cell_pts_ids = Tree.cell_pts(parapet_cell_id).id;
            parapet_cell_pts_xyz = Tree.pts(parapet_cell_pts_ids,1:3);
            
            % Search cells on the footpath
            mask = ismember(footpath_cell_ids(:,1), parapet_cell_id);
            if any(mask)
                % Available cells in the footpath, which is used as a local
                % surface of the footpath
                footpath_parapet_cell_ptc_ids = supstruct(1).Component(footpath_id).cell(mask).ptc_ids;
                footpath_parapet_cell_surf = supstruct(1).Component(footpath_id).cell(mask).surface_features;

            else
                % No available so footpath cells adjacent to the cell will
                % be used to determine a local surface of the footpath 
                flag = true;
                scale = 1.0;
                while flag
                    neighbour_cell_ids = Window_Query_Neighbour_Cell(Tree, parapet_cell_id, footpath_cell_ids(:,1), scale*threshold.cell_size);
                    if isempty(neighbour_cell_ids)
                        scale = scale + 1;
                    else
                        flag = false;
                    end   
                end
                mask = ismember(footpath_cell_ids(:,1),neighbour_cell_ids);
                footpath_parapet_cell_ptc_ids = vertcat(supstruct(1).Component(footpath_id).cell(mask).ptc_ids);
                footpath_parapet_cell_ptc_xyz = Tree.pts(footpath_parapet_cell_ptc_ids,1:3);
                footpath_parapet_cell_ptc_cent = mean(footpath_parapet_cell_ptc_xyz,1);
                footpath_parapet_cell_ptc_normal = eigenspace(footpath_parapet_cell_ptc_xyz,1);
                footpath_parapet_cell_surf = [footpath_parapet_cell_ptc_cent, footpath_parapet_cell_ptc_normal];
                clear flag scale footpath_parapet_cell_ptc_xyz footpath_parapet_cell_ptc_cent footpath_parapet_cell_ptc_normal
            end
            
            % Extract the point above the footpath
            if dot(footpath_parapet_cell_surf(4:6), threshold.nz) < 0
                footpath_parapet_cell_surf(4:6) = - footpath_parapet_cell_surf(4:6);
            end
            dist_pts_footpath = dist_3Dpoints_3Dplane(parapet_cell_pts_xyz, footpath_parapet_cell_surf);
            mask = 2.0*threshold.max_distance <= dist_pts_footpath;
            dist_pts_footpath = dist_pts_footpath(mask);
            parapet_cell_pts_ids = parapet_cell_pts_ids(mask);
            parapet_cell_pts_xyz = parapet_cell_pts_xyz(mask,:);
            
            % Remove the points assigned to the footpath if they are available
            if numel(parapet_cell_pts_ids) <= threshold.min_num_pts
                continue;
            end
            mask = ismember(parapet_cell_pts_ids,footpath_parapet_cell_ptc_ids);
            dist_pts_footpath(mask) = [];
            parapet_cell_pts_ids(mask) = [];
            parapet_cell_pts_xyz(mask,:) = [];
            
            % Update the structure
            if (numel(parapet_cell_pts_ids) >= threshold.min_num_pts)&&(max(dist_pts_footpath) >= struct_threshold.parapet_min_height)
                % Update parapet_cell
                parapet.cell(no_parapet_cell).id = [parapet_cell_id, 0];
%                 cell_mbb = min2DBoundingBox(parapet_cell_ptc_xyz(:,1:2)');
%                 parapet.cell(no_parapet_cell).bounds = mean(cell_mbb.vertices,1);
                parapet.cell(no_parapet_cell).bounds = [min(parapet_cell_pts_xyz,[],1),max(parapet_cell_pts_xyz,[],1)];
                parapet.cell(no_parapet_cell).ptc_ids = parapet_cell_pts_ids;
                parapet.cell(no_parapet_cell).surface_features = [mean(parapet_cell_pts_xyz,1), eigenspace(parapet_cell_pts_xyz,1)];
                parapet.cell(no_parapet_cell).height =  max(dist_pts_footpath);
                no_parapet_cell = no_parapet_cell + 1;
            end
            clear parapet_cell_pts_ids parapet_cell_pts_xyz dist_pts_footpath footpath_parapet_cell_surf
        
        end
        
        %% Classify vehicle parapet and pedestrian parapet       
        % Preallocation
        final_parapet = struct('cell',[],'status',[]);
        num_parapet = 1;
        % Retrieve features of the cells
        parapet_cell_ids = vertcat(parapet.cell.id);
        parapet_cell_height = vertcat(parapet.cell.height);
%         parapet_cell_bounds = vertcat(parapet.cell.bounds);
        parapet_cell_bounds_cent = vertcat(parapet.cell.surface_features);
        parapet_cell_bounds_cent = parapet_cell_bounds_cent(:,1:3);
%         parapet_cell_bounds_cent = 0.5*(parapet_cell_bounds(:,[1,2,3]) + parapet_cell_bounds(:,[4,5,6]));
        parapet_cell_bounds_proj_cent = proj_3Dpoints_3Dplane(parapet_cell_bounds_cent, [road_surface.cent, road_surface.normal]);
        clear parapet_cell_bounds parapet_cell_bounds_cent 
        
        % Calculate distance from the cell to the road central line
%         dist_parapet_road = abs(dist_3Dpoints_3Dline(parapet_cell_bounds_proj_cent, [road_surface.cent,road_surface.tangent]));
        dist_parapet_road = abs(dist_3Dpoints_3Dline(parapet_cell_bounds_proj_cent, [road_surface.cent,road_surface.mbb.long_edge_vector]));
        
        % Use a histogram based distance to classify the cells: vehicle and pedestrian parapets
        num_bin = ceil((max(dist_parapet_road) - min(dist_parapet_road))/bin_width_dist);
        dist_bin_edges = 0.5*(max(dist_parapet_road) + min(dist_parapet_road) - num_bin*bin_width_dist):bin_width_dist:0.5*(max(dist_parapet_road) + min(dist_parapet_road) + num_bin*bin_width_dist);
        [~, ~, dist_cell_bin_ids] = histcounts(dist_parapet_road, dist_bin_edges);
        clear num_bin dist_bin_edges
        
        % Filtering the short bin based on a projected length of the cells
        % within the bin on the road central line
        dist_bin_ids = unique(dist_cell_bin_ids);
        dist_bin_info = inf(numel(dist_bin_ids),2);
        dist_bin_info(:,1) = dist_bin_ids; 
        for k = 1:numel(dist_bin_ids)
            % Retrieve the cells within the bin
            dist_bin_id = dist_bin_ids(k);
            mask = ismember(dist_cell_bin_ids,dist_bin_id);
            dist_cell_road = dist_parapet_road(mask);
            dist_cell_proj_cent = parapet_cell_bounds_proj_cent(mask,:);
            
            % Cal projected length
            [~, dist_proj_length] = cal_proj_length_3Dpoints(dist_cell_proj_cent, [road_surface.cent,road_surface.mbb.long_edge_vector], []);
            if (struct_threshold.bridge_min_length <= dist_proj_length)&&(min_num_cells <= numel(dist_cell_proj_cent))
                dist_bin_info(k,2) = mean(dist_cell_road);
            end
            clear dist_bin_id dist_cell_road dist_cell_proj_cent dist_proj_length
        end
        mask = any(isinf(dist_bin_info),2);
        dist_bin_info = dist_bin_info(~mask, :);

        % Cluster the bins based on a distance between them
        diff_dist_cell_road = diff(dist_bin_info(:,2));
        mask = 0.75*struct_threshold.footpath_min_width < diff_dist_cell_road;
        local_ids = find(mask == 1);
        dist_bin_local_ids = [0, local_ids, size(dist_bin_info,1)];
        clear diff_dist_cell_road local_ids
        
        % Group the bins
        cluster_dist_bin = zeros(numel(dist_bin_local_ids) -1, 3);
        cluster_dist_bin(:,1) = 1:numel(dist_bin_local_ids) - 1;
        cluster_dist_bin(:,2) = dist_bin_info(dist_bin_local_ids(1:end-1)+1,1);
        cluster_dist_bin(:,3) = dist_bin_info(dist_bin_local_ids(2:end),1);
        
        % Use the histogram based on the height 
        for k = 1:size(cluster_dist_bin,1)
            % Retrieve the cells in the cluster
            cluster_ids = cluster_dist_bin(k,2):1:cluster_dist_bin(k,3);
            mask = ismember(dist_bin_info(:,1), cluster_ids);
            cluster_dist_bin_ids = dist_bin_info(mask,1);
            mask = ismember(dist_cell_bin_ids, cluster_dist_bin_ids);
            cluster_dist_bin_cell_ids = parapet_cell_ids(mask, :);
            cluster_dist_bin_cell_height = parapet_cell_height(mask);
            
            % Use dbscan to remove outlier cell
            db_cluster = dbscan(cluster_dist_bin_cell_height,struct_threshold.parapet_height_tol,ceil(0.5*min_num_cells));
            mask = db_cluster > 0;
            filter_cell_ids = cluster_dist_bin_cell_ids(mask, :);
            
%             mask = db_cluster <= 0;
%             cluster_dist_bin_cell_height(mask)
%             cluster_dist_bin_cell_ids(mask)
%             
%             cluster_dist_bin_cell_height(db_cluster == 2)

%             % Use a histogram based cells' height
%             bin_width_height = 0.5*mean(cluster_dist_bin_cell_height);
%             num_bin = ceil((max(cluster_dist_bin_cell_height) - min(cluster_dist_bin_cell_height))/bin_width_height);
%             height_bin_edges = 0.5*(max(cluster_dist_bin_cell_height) + min(cluster_dist_bin_cell_height) - num_bin*bin_width_height):bin_width_height:0.5*(max(cluster_dist_bin_cell_height) + min(cluster_dist_bin_cell_height) + num_bin*bin_width_height);
%             [height_bin_count, ~, height_bin_cell_ids] = histcounts(cluster_dist_bin_cell_height, height_bin_edges);
%             clear bin_width_height num_bin height_bin_edges 
%             
%             % Filtering the outlier cells
%             mask = 0.25*numel(cluster_dist_bin_cell_height) <= height_bin_count;
%             height_bin_local_ids = find(mask);
%             mask = ismember(height_bin_cell_ids, height_bin_local_ids);
%             filter_cell_ids = cluster_dist_bin_cell_ids(mask, :); 
%             clear height_bin_local_ids
            
            % Update the data structure
            mask = ismember(parapet_cell_ids, filter_cell_ids, 'rows');
            cell_local_ids = find(mask == 1);
            
            final_parapet(num_parapet).status = 1;
            for count = 1:numel(cell_local_ids)
                local_id = cell_local_ids(count);
                final_parapet(num_parapet).cell(count).id = parapet.cell(local_id).id;
                final_parapet(num_parapet).cell(count).ptc_ids = parapet.cell(local_id).ptc_ids;
                final_parapet(num_parapet).cell(count).surface_features = parapet.cell(local_id).surface_features;
            end
            num_parapet = num_parapet + 1;
        end
        
        %% Update the bridge tree
        if length(final_parapet) == 1
            comp_name = {'Vehicle Parapet'};
        elseif length(final_parapet) == 2 
            comp_name = {{'Vehicle Parapet'}, {'Pedestrian Parapet'}};
        else
            comp_name = cell(1,length(final_parapet));
            comp_name(:) = {'Parapet'};
        end
        for count = 1:length(final_parapet)
            % Retrieve name of the footpath
            footpath_name = supstruct(1).Description{footpath_id};
            footpath_num = sscanf(footpath_name,'Footpath%d');
            slave_comp_id = count;
            slave_comp_name = strcat(comp_name{count}, " ", num2str(footpath_num));
            [supstruct, ~] = bridge_tree_components(supstruct, footpath_name, slave_comp_id, slave_comp_name, final_parapet);
        end
    end
end
end % End of the function

%%
% Preallocation
%         final_parapet = struct('cell',[],'status',[]);
%         num_parapet = 1;
%         % Retrieve features of the cells
%         parapet_cell_ids = vertcat(parapet.cell.id);
%         parapet_cell_height = vertcat(parapet.cell.height);
%         parapet_cell_bounds = vertcat(parapet.cell.bounds);
%         parapet_cell_bounds_cent = 0.5*(parapet_cell_bounds(:,[1,2,3]) + parapet_cell_bounds(:,[4,5,6]));
%         parapet_cell_bounds_proj_cent = proj_3Dpoints_3Dplane(parapet_cell_bounds_cent, [road_surface.cent, road_surface.normal]);
%         clear parapet_cell_bounds parapet_cell_bounds_cent 
%         
%         % Calculate distance from the cell to the road central line
%         dist_parapet_road = abs(dist_3Dpoints_3Dline(parapet_cell_bounds_proj_cent, [road_surface.cent,road_surface.tangent]));
%         
%         % Use a histogram based distance to classify the cells: vehicle and
%         % pedestrian parapets
%         num_bin = ceil((max(dist_parapet_road) - min(dist_parapet_road))/bin_width_dist);
%         dist_bin_edges = 0.5*(max(dist_parapet_road) + min(dist_parapet_road) - num_bin*bin_width_dist):bin_width_dist:0.5*(max(dist_parapet_road) + min(dist_parapet_road) + num_bin*bin_width_dist);
%         [~, ~, dist_cell_bin_ids] = histcounts(dist_parapet_road, dist_bin_edges);
%         clear num_bin dist_bin_edges
%         
%         % Filtering the short bin based on a projected length of the cells
%         % within the bin on the road central line
%         dist_bin_ids = unique(dist_cell_bin_ids);
%         dist_bin_info = inf(numel(dist_bin_ids),2);
%         dist_bin_info(:,1) = dist_bin_ids; 
%         for k = 1:numel(dist_bin_ids)
%             % Retrieve the cells within the bin
%             dist_bin_id = dist_bin_ids(k);
%             mask = ismember(dist_cell_bin_ids,dist_bin_id);
%             dist_cell_road = dist_parapet_road(mask);
%             dist_cell_proj_cent = parapet_cell_bounds_proj_cent(mask,:);
%             
%             % Cal projected length
%             [~, dist_proj_length] = cal_proj_length_3Dpoints(dist_cell_proj_cent, [road_surface.cent,road_surface.tangent], []);
%             if (struct_threshold.bridge_min_length <= dist_proj_length)&&(min_num_cells <= numel(dist_cell_proj_cent))
%                 dist_bin_info(k,2) = mean(dist_cell_road);
%             end
%             clear dist_bin_id dist_cell_road dist_cell_proj_cent dist_proj_length
%         end
%         mask = any(isinf(dist_bin_info),2);
%         dist_bin_info = dist_bin_info(~mask, :);
% %         [~, sort_ids] = sort(dist_bin_info(:,2)); % Data sort already
% %         dist_bin_info = dist_bin_info(sort_ids,:);
%         
%         % Cluster the bins based on a distance between them
%         diff_dist_cell_road = diff(dist_bin_info(:,2));
%         mask = 0.75*struct_threshold.footpath_min_width < diff_dist_cell_road;
%         local_ids = find(mask == 1);
%         dist_bin_local_ids = [0, local_ids, size(dist_bin_info,1)];
%         clear diff_dist_cell_road local_ids
%         
%         % Group the bins
%         cluster_dist_bin = zeros(numel(dist_bin_local_ids) -1, 3);
%         cluster_dist_bin(:,1) = 1:numel(dist_bin_local_ids) - 1;
%         cluster_dist_bin(:,2) = dist_bin_info(dist_bin_local_ids(1:end-1)+1,1);
%         cluster_dist_bin(:,3) = dist_bin_info(dist_bin_local_ids(2:end),1);
%         
%         % Use the histogram based on the height 
%         for k = 1:size(cluster_dist_bin,1)
%             % Retrieve the cells in the cluster
%             cluster_ids = cluster_dist_bin(k,2):1:cluster_dist_bin(k,3);
%             mask = ismember(dist_bin_info(:,1), cluster_ids);
%             cluster_dist_bin_ids = dist_bin_info(mask,1);
%             mask = ismember(dist_cell_bin_ids, cluster_dist_bin_ids);
%             cluster_dist_bin_cell_ids = parapet_cell_ids(mask, :);
%             cluster_dist_bin_cell_height = parapet_cell_height(mask);
% 
%             % Use a histogram based cells' height
%             bin_width_height = 0.5*mean(cluster_dist_bin_cell_height);
%             num_bin = ceil((max(cluster_dist_bin_cell_height) - min(cluster_dist_bin_cell_height))/bin_width_height);
%             height_bin_edges = 0.5*(max(cluster_dist_bin_cell_height) + min(cluster_dist_bin_cell_height) - num_bin*bin_width_height):bin_width_height:0.5*(max(cluster_dist_bin_cell_height) + min(cluster_dist_bin_cell_height) + num_bin*bin_width_height);
%             [height_bin_count, ~, height_bin_cell_ids] = histcounts(cluster_dist_bin_cell_height, height_bin_edges);
%             clear bin_width_height num_bin height_bin_edges 
%             
%             % Filtering the outlier cells
%             mask = 0.25*numel(cluster_dist_bin_cell_height) <= height_bin_count;
%             height_bin_local_ids = find(mask);
%             mask = ismember(height_bin_cell_ids, height_bin_local_ids);
%             filter_cell_ids = cluster_dist_bin_cell_ids(mask, :); 
%             clear height_bin_local_ids
%             
%             % Update the data structure
%             mask = ismember(parapet_cell_ids, filter_cell_ids, 'rows');
%             cell_local_ids = find(mask == 1);
%             
%             final_parapet(num_parapet).status = 1;
%             for count = 1:numel(cell_local_ids)
%                 local_id = cell_local_ids(count);
%                 final_parapet(num_parapet).cell(count).id = parapet.cell(local_id).id;
%                 final_parapet(num_parapet).cell(count).ptc_ids = parapet.cell(local_id).ptc_ids;
%                 final_parapet(num_parapet).cell(count).surface_features = parapet.cell(local_id).surface_features;
%             end
%             num_parapet = num_parapet + 1;
%         end
% Retrieve features of the cells
% parapet_cell_ids = vertcat(parapet.cell.id);
% parapet_cell_height = vertcat(parapet.cell.height);
% parapet_cell_bounds = vertcat(parapet.cell.bounds);
% parapet_cell_bounds_cent = 0.5*(parapet_cell_bounds(:,[1,2,3]) + parapet_cell_bounds(:,[4,5,6]));
% parapet_cell_bounds_proj_cent = proj_3Dpoints_3Dplane(parapet_cell_bounds_cent, [road_surface_cent, road_surface_normal]);
% dist_parapet_road = abs(dist_3Dpoints_3Dline(parapet_cell_bounds_proj_cent, [road_surface_cent,road_surface_tangent]));
% 
% % Use distance
% [fi,zi,~] = ksdensity(dist_parapet_road,'npoints',100,'bandwidth',0.2*struct_threshold.footpath_min_width,'Kernel','epanechnikov');
% peak_shape_dist = peak_shape_width(zi, fi);
% 
% figure(1)
% %         plot(zi, fi)
% h = histogram(dist_parapet_road)
% 
% 
% % Preallocation
% final_parapet = struct('cell',[],'status',[]);
% num_parapet = 1;
% for k = 1:size(peak_shape_dist,1)
%     % Check peak k
%     mask = (peak_shape_dist(k,1) - peak_shape_dist(k,2) <= dist_parapet_road)&(dist_parapet_road <= peak_shape_dist(k,1) + peak_shape_dist(k,3));
%     dist_peak_parapet_cell_height = parapet_cell_height(mask,:);
%     dist_peak_parapet_cell_ids = parapet_cell_ids(mask,:);
% 
%     % use height
%     [fi,zi,~] = ksdensity(dist_peak_parapet_cell_height(:,2),'npoints',100,'bandwidth',0.5*struct_threshold.parapet_min_height,'Kernel','epanechnikov');
%     peak_shape_dist_height = peak_shape_width(zi, fi);
% 
%     for kk = 1:size(peak_shape_dist_height,1)
%         mask = (peak_shape_dist_height(kk,1) - peak_shape_dist_height(kk,2) <= dist_peak_parapet_cell_height(:,2))&...
%                (dist_peak_parapet_cell_height(:,2) <= peak_shape_dist_height(kk,1) + peak_shape_dist_height(kk,3));
% 
%         dist_height_peak_parapet_cell_ids = dist_peak_parapet_cell_ids(mask,:);
% 
%         % Compute the distance of the cluster
%         mask = ismember(parapet_cell_ids, dist_height_peak_parapet_cell_ids, 'row');
%         cell_local_ids = find(mask ==1);
%         dist_height_peak_parapet_cell_proj_cent = parapet_cell_bounds_proj_cent(cell_local_ids,:);
%         [~, proj_length] = cal_proj_length_3Dpoints(dist_height_peak_parapet_cell_proj_cent, [road_surface_cent,road_surface_tangent], []);
% 
%         % Update the final structure
%         if (proj_length >= struct_threshold.bridge_min_length)&&(struct_threshold.bridge_min_length/threshold.cell_size <= numel(cell_local_ids))
% 
%             final_parapet(num_parapet).status = 1;
%             for count = 1:numel(cell_local_ids)
%                 local_id = cell_local_ids(count);
%                 final_parapet(num_parapet).cell(count).id = parapet.cell(local_id).id;
%                 final_parapet(num_parapet).cell(count).ptc_ids = parapet.cell(local_id).ptc_ids;
%             end
%             num_parapet = num_parapet + 1;
% 
%         end
%     end
% 
% 
% 
% end