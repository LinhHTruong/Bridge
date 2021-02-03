function [supstruct, active_pts] = bridge_intermediate_superstructure(Tree, supstruct, struct_threshold, threshold, time_record)
%% This function is to extract the components surface of the bridge structure
% Input:
%     Tree                      : Data structure stored 2D cell in xy associated with ptc_ids
%     super_structure           : Data structure store results
%     struct_threshold          : Structur_threshold relating a minimum size of the structural components
%     threshold                 : Threashold relating to filtering and rocessing point cloud
%     time_record               : Check running time

% Output
%     super_structure           : Data structure store results
%     active_pts                : [Nx1] - True and False marks the points still active

%% Demo
% Tree = OQTR;
% superstruct = SuperStructure;
% struct_threshold;
% threshold = THRESHOLD;
% time_record = true;
if time_record
    tic
end
%% Set up link_list for the points within
active_pts = true(size(Tree.pts,1),1);
%% Extract top surfaces of a supperstructure
% Extract the road surface
mask = cellfun(@(s) contains(s, 'Road Surface', 'IgnoreCase', true), supstruct(1).Description, 'UniformOutput', false);
road_surf_id = find(cell2mat(mask));

% Extract footpath
mask = cellfun(@(s) contains(s, 'Footpath', 'IgnoreCase', true), supstruct(1).Description, 'UniformOutput', false);
footpath_ids = find(cell2mat(mask));

% Merge top surfaces of a superstructure
supstruct_top_surf = supstruct.Component(road_surf_id).cell;
for i = 1:numel(footpath_ids)
    supstruct_top_surf = [supstruct_top_surf,supstruct.Component(footpath_ids(i)).cell]; 
end

% Extract beam bottom
mask = cellfun(@(s) strcmp(s, 'Bottom Beam'), supstruct(1).Description, 'UniformOutput', false);
beam_main_bot_surf_id = find(cell2mat(mask));

mask = cellfun(@(s) contains(s, 'Beam', 'IgnoreCase', true), supstruct(1).Description, 'UniformOutput', false);
beam_bot_surf_ids = find(cell2mat(mask));
beam_bot_surf_ids = setdiff(beam_bot_surf_ids,beam_main_bot_surf_id); 
% Merge top surfaces of a superstructure
supstruct_bot_surf = supstruct.Component(beam_main_bot_surf_id).cell;
for i = 1:numel(beam_bot_surf_ids)
    supstruct_bot_surf = [supstruct_bot_surf,supstruct.Component(beam_bot_surf_ids(i)).cell]; 
end

%% Extract the cells within the 
intermediate_surf = struct('cell',[]);
% Retrieve cell_peaks on the top surface
supstruct_top_cell_peak_ids = vertcat(supstruct_top_surf.ids);
supstruct_top_cell_ids = unique(supstruct_top_cell_peak_ids(:,1));
supstruct_top_cell_surf = vertcat(supstruct_top_surf.surface_features);

% Retrieve cell_peaks on the bottom surface
supstruct_bot_cell_peak_ids = vertcat(supstruct_bot_surf.ids);
supstruct_bot_cell_ids = unique(supstruct_bot_cell_peak_ids(:,1));
supstruct_bot_cell_surf = vertcat(supstruct_bot_surf.surface_features);
no_cells = 1;

% Cells of superstructure
supstruct_bot_cell_ids = union(supstruct_top_cell_ids,supstruct_bot_cell_ids); 
for i = 1:numel(supstruct_bot_cell_ids)
    
    % Retrieve the cell_peaks on the top and bottom surface
    cell_id = supstruct_bot_cell_ids(i);
    top_mask = ismember(supstruct_top_cell_peak_ids(:,1), cell_id);
    bot_mask = ismember(supstruct_bot_cell_peak_ids(:,1), cell_id);
    
    if any(top_mask)&& any(bot_mask) 
        % Get the top surface directly
        [top_surf_cell_ptc_ids, top_surf_cell_surf] = get_cell_surf(supstruct_top_surf, top_mask, 'top');
        
        % Get the bottom surface directly
        [bot_surf_cell_ptc_ids, bot_surf_cell_surf] = get_cell_surf(supstruct_bot_surf, bot_mask, 'bottom');
        
    elseif ~any(top_mask) && any(bot_mask)
        % Get the bottom surface directly
        [bot_surf_cell_ptc_ids, bot_surf_cell_surf] = get_cell_surf(supstruct_bot_surf, bot_mask, 'bottom'); 
        
        % Searching the neighbour cells of the bottom surface in the top surfaces for the cell 
        neighbour_cell_id = knnsearch(supstruct_top_cell_surf(:,[1,2]), bot_surf_cell_surf(:,[1,2]),'K',1);
        top_surf_cell_ptc_ids = [];
        top_surf_cell_surf = supstruct_top_cell_surf(neighbour_cell_id,:);
    else %any(top_mask) && ~any(bot_mask)
        
        % Get the top surface directly
        [top_surf_cell_ptc_ids, top_surf_cell_surf] = get_cell_surf(supstruct_top_surf, top_mask, 'top');
        
        % Searching the neighbour cells of the top surface in the bottom surfaces for the cell 
        neighbour_cell_id = knnsearch(supstruct_bot_cell_surf(:,[1,2]), top_surf_cell_surf(:,[1,2]),'K',1);
        bot_surf_cell_ptc_ids = [];
        bot_surf_cell_surf = supstruct_bot_cell_surf(neighbour_cell_id,:);
        
    end
    
    % Adjust direction of the surface
    if dot(top_surf_cell_surf(4:6), threshold.nz) < 0
        top_surf_cell_surf(4:6) = -top_surf_cell_surf(4:6);
    end
    if dot(bot_surf_cell_surf(4:6), threshold.nz) < 0
        bot_surf_cell_surf(4:6) = -bot_surf_cell_surf(4:6);
    end
    
    % Retrieve the points within the cell
    cell_pts_ids = Tree.cell_pts(cell_id).id;
    
    % Deactive all points in the cell
    active_pts(cell_pts_ids) = false;
    
    % Remove the points were assigned to the top and bottom
    cell_pts_ids = setdiff(cell_pts_ids,top_surf_cell_ptc_ids); 
    cell_pts_ids = setdiff(cell_pts_ids,bot_surf_cell_ptc_ids);
    
    % Retrieve the remaining points
    cell_pts_xyz = Tree.pts(cell_pts_ids, 1:3);
    
    % Extract the points between
    % Cal distance from the points to the surfac
    dist_ptc_top_surface = dist_3Dpoints_3Dplane(cell_pts_xyz, top_surf_cell_surf);
    dist_ptc_bot_surface = dist_3Dpoints_3Dplane(cell_pts_xyz, bot_surf_cell_surf);
    
    % Retrieve the intermediate points
    mask = (dist_ptc_top_surface <= -threshold.max_residual)&(threshold.max_residual <= dist_ptc_bot_surface);
    intermediate_cell_pts_ids = cell_pts_ids(mask);
    
    if numel(intermediate_cell_pts_ids) >= threshold.min_num_pts
        % Retrieve the points coordinates    
        intermediate_cell_pts_xyz = cell_pts_xyz(mask,:);
    
        % Compute the features
        intermediate_cell_cent = mean(intermediate_cell_pts_xyz,1);
        intermediate_cell_norm = eigenspace(intermediate_cell_pts_xyz,1);

        % Compute the 
        intermediate_surf.cell(no_cells).id = [cell_id, 0];
        intermediate_surf.cell(no_cells).ptc_ids = intermediate_cell_pts_ids;
        intermediate_surf.cell(no_cells).surface_features = [intermediate_cell_cent, intermediate_cell_norm];
        no_cells = no_cells + 1;
    end
    
    % Active all the points below the bottom surface
    mask = dist_ptc_bot_surface <= 0;
    cell_lower_ptc_ids = cell_pts_ids(mask);
    active_pts(cell_lower_ptc_ids) = true;
end

%% Determine traffice direction
road_surf_pts_ids = vertcat(supstruct.Component(road_surf_id).cell.ptc_ids);
road_surf_pts_xyz = Tree.pts(road_surf_pts_ids,1:3);
road_mbb = min2DBoundingBox(road_surf_pts_xyz(:,1:2)');
road_surf_cent = [mean(road_mbb.vertices,1), 0];
road_surf_normal = cross([road_mbb.long_edge_vector, 0], [road_mbb.short_edge_vector,0]);
road_surf_tangent = [road_mbb.long_edge_vector,0];

%% Filtering unrealistic cells by using cell connectivity
% Retrieve the cell center
cells_cent = vertcat(intermediate_surf.cell.surface_features);
cells_cent = cells_cent(:,1:3);

% Compute the distance from the cells to the traffic direction
cells_proj_cent = proj_3Dpoints_3Dplane(cells_cent, [road_surf_cent, road_surf_normal]);
dist_cells_road = dist_3Dpoints_3Dline(cells_proj_cent, [road_surf_cent,road_surf_tangent], cross(road_surf_normal, road_surf_tangent));
        
%% Use a histogram based distance to classify the cells: vehicle and
% pedestrian parapets
[fi,zi,~] = ksdensity(dist_cells_road,'npoints',100,'bandwidth',0.05,'Kernel','epanechnikov');
peak_shape_dist = peak_shape_width(zi, fi);
%% Cluster the cell within a peak
% intermediate_ptc_ids = vertcat(intermediate_surf.cell.ptc_ids);
% flag = false(numel(intermediate_ptc_ids),1);
intermediate_seg = struct('ptc_ids',[], 'cell_ids',[]);
no_seg = 1;
voxel_region_growing_threshold = voxel_region_growing_params(1, 0.2, threshold.max_angle, threshold.sampling_step);
for i = 1:size(peak_shape_dist)
    % Retrieve the cells within the peak
    mask = (peak_shape_dist(i,1) - peak_shape_dist(i,2) <= dist_cells_road) & (dist_cells_road <= peak_shape_dist(i,1) + peak_shape_dist(i,3));
    cell_peak_ids = vertcat(intermediate_surf.cell(mask).id);
    cell_peak_pts_ids = vertcat(intermediate_surf.cell(mask).ptc_ids);
    cell_peak_pts_xyz = Tree.pts(cell_peak_pts_ids,1:3);
    
    % Voxel-based segmentation
    ptc_segment_info = voxel_region_growing_segmentation(cell_peak_pts_xyz,voxel_region_growing_threshold);
    
    % Merging the planar segments
    ptc_segment_info = plane_segment_merging(cell_peak_pts_xyz, ptc_segment_info, struct_threshold, threshold, false);
    
    % Eliminate the small segment
    ptc_segment_info = planar_surface_filtering(cell_peak_pts_xyz, ptc_segment_info, 'length', struct_threshold, threshold, false);

    % Update the inlier points
    if ~isempty(ptc_segment_info)
        intermediate_seg(no_seg).ptc_ids = cell_peak_pts_ids(ptc_segment_info(:,1));
        intermediate_seg(no_seg).cell_ids = cell_peak_ids;
        no_seg = no_seg + 1;
    end
    clear cell_peak_pts_ids cell_peak_pts_xyz ptc_segment_info
end

%% Update the cells
% Retrieve the cell ids of the intermediate 
cell_ids = vertcat(intermediate_surf.cell.id);

% Preallocation 
final_intermediate_surf = struct('cell',[], 'status', []);

% Update process
for i = 1:length(intermediate_seg)
    % Retrieve the points within the segment
    seg_cell_pts_ids = intermediate_seg(i).ptc_ids;
    seg_cell_ids = intermediate_seg(i).cell_ids;
    
    % Retrieve the cell
    no_cell = 1;
    for j = 1:size(seg_cell_ids, 1)
        
        % Retrieve the points within the cell
        mask = ismember(cell_ids(:,1),seg_cell_ids(j,1));
        cell_pts_ids = intermediate_surf.cell(mask).ptc_ids;
        cell_surf_features = intermediate_surf.cell(mask).surface_features;

        % Check if points within the cells are inlier
        mask = ismember(cell_pts_ids,seg_cell_pts_ids);
    
        if all(mask)
            % All points within the cell are inlier
            final_intermediate_surf(i).cell(no_cell).id = seg_cell_ids(j,:);
            final_intermediate_surf(i).cell(no_cell).ptc_ids = cell_pts_ids;
            final_intermediate_surf(i).cell(no_cell).surface_features = cell_surf_features;
             no_cell = no_cell + 1;
        else
            if sum(mask) > threshold.min_num_pts
                % Parts of the points within the cell are inliers
                cell_pts_xyz = Tree.pts(cell_pts_ids,1:3);
                final_intermediate_surf(i).cell(no_cell).id = seg_cell_ids(j,:);
                final_intermediate_surf(i).cell(no_cell).ptc_ids = cell_pts_ids;
                final_intermediate_surf(i).cell(no_cell).surface_features = [mean(cell_pts_xyz,1), eigenspace(cell_pts_xyz,1)];
                no_cell = no_cell + 1;
            end
        end
    end
    
    % Update status
    final_intermediate_surf(i).status = true;
end

%% Update superstructure
for i = 1:length(final_intermediate_surf)
    
    % Master segment
    slave_comp_ids = i;
    slave_comp_name = strcat('Intermediate Surface', " ", num2str(i));
    
    % Assign the data structure
    [supstruct, ~] = bridge_tree_components(supstruct, ' ', slave_comp_ids, slave_comp_name, final_intermediate_surf);
    
end
if time_record
    fprintf('Running time of extracting intermediate surfaces of a superstructure: %.2f seconds \n', toc);
end