function [supstruct, cells_region, un_connect_segs_ids] = bridge_bottom_superstructure(Tree, cells_region, supstruct, struct_threshold, threshold, time_record)
%% This function is to extract the components surface of the bridge structure
% Input:
%     Tree                      : Data structure stored 2D cell in xy associated with ptc_ids
%     bridge_tree               : Data structure store results
%     region                    : Data structure stored points in each segment (plane)
%     comp_name                 : Name of supperstructure
%     struct_threshold          : Structur_threshold relating a minimum size of the structural components
%     threshold                 : Threashold relating to filtering and rocessing point cloud
%     debug                     : Check running time

% Output
%     bridge_tree               : Data structure store results
%     region                    : Data structure stored points in each segment (plane)
%     comp_link_list            : [Nx4] [ref_seg, connect_seg, relationship, side]; 
%                                 relationship: coplanar, parallel, intersection; side = 1 and -1

%% Demo
%%
% Tree = OQTR;
% supstruct = SuperStructure;
% cells_region = Cells_Region;
% threshold = THRESHOLD;
% time_record = true;

% Return results
% SuperStructure = supstruct;
% Cells_Region = cells_region;
% un_connect_segs_ids = un_connect_segs_ids;
%% Compute the feature of the region
region_features = cal_region_features(Tree, cells_region);
%% Find the reference component: bottom surface of girder
% Retrieve the cells and its surfaces of the road surface
mask = cellfun(@(s) contains('Road Surface', s, 'IgnoreCase', true), supstruct(1).Description);
% road_surf_id = find(mask);
road_surf_cell_ids = vertcat(supstruct.Component(mask).cell.ids);
road_surf_cell_surf = vertcat(supstruct.Component(mask).cell.surface_features);

% Compute distance from the bottom surfaces of the superstructure to
% the road surface by using cell-to-cell method
no_region = length(cells_region);
dist_region_road_surface = [];
for region_id = 1:no_region
    % Retrieve the cells within the region
    cell_ids = vertcat(cells_region(region_id).cell.id); 
    cell_surf_features = vertcat(cells_region(region_id).cell.surface_features); 

    % Compute a distance from region to road surface
    dist_cell_road_surface = zeros(size(cell_surf_features,1), 1);

    for i = 1:size(cell_ids,1)
        % Search overlap cells
        mask = ismember(road_surf_cell_ids(:,1), cell_ids(i,1));
        if any(mask)
            % Retrieve a local surface on the bridge road surface: the
            % same cell id
            bridge_road_local_surf = road_surf_cell_surf(mask,1:6);
        else
            % Retrieve a local surface on the bridge road surface: the
            % neastest cells
            neighbour_cell_id = knnsearch(road_surf_cell_surf(:,[1,2]),cell_surf_features(i,[1,2]),'K',1);
            bridge_road_local_surf = road_surf_cell_surf(neighbour_cell_id,1:6);
            clear neighbour_cell_id
        end

        % Compute distance: delta_z
%             intersection_pt = cal_line_plane_intersection([cell_surf_features(i,1:3),[0,0,1]], bridge_road_local_surf);
%             proj_ptc = proj_3Dptc_3Dplane(cell_surf_features(i,1:3),bridge_road_local_surf)
        dist_cell_road_surface(i) = abs(cell_surf_features(i,3) - bridge_road_local_surf(3));%abs(dist_3Dpoints_3Dplane(cell_surf_features(i,1:3),bridge_road_local_surf));
        clear mask bridge_road_local_surf
    end
    % Assembly data
%         region_cell_dist = [cell_surf_features(:,[1,2,3]),dist_cell_road_surface, cell_ids];
    dist_cell_road_surface(:,2) = region_id;
    dist_region_road_surface = [dist_region_road_surface;dist_cell_road_surface];
    clear cell_ids cell_surf_features dist_cell_road_surface
end 
clear bridge_comp_name bridge_road_surf_id bridge_road_surf_cell_ids bridge_road_surf_cell_surf
%% Determine the regions are the bottom surface of the bridge desk by using kde: the distance is mostly constant: DBscan is another option for curvature
[fi,zi,~] = ksdensity(dist_region_road_surface(:,1),'npoints',100,'bandwidth',threshold.min_dist_2_planes,'Kernel','epanechnikov');
peak_shape = peak_shape_width(zi, fi);
% Count the number of the cell in the peaks
peak_count_cells = zeros(size(peak_shape,1),1);
for count = 1:size(peak_shape,1)
    mask = (peak_shape(count,1) - peak_shape(count,2) <= dist_region_road_surface(:,1))&...
           (dist_region_road_surface(:,1) <= peak_shape(count,1) + peak_shape(count,3)); 
    peak_count_cells(count) = sum(mask);
end

[~, max_id] = max(peak_count_cells);
mask = (peak_shape(count,1) - peak_shape(max_id,2) <= dist_region_road_surface(:,1))&...
       (dist_region_road_surface(:,1) <= peak_shape(count,1) + peak_shape(max_id,3));   
largest_region_cell_ids = dist_region_road_surface(mask,2);
clear fi zi peak_shape max_id mask peak_count_cells
%% Decision which regions are the bottom of the bridge desk
% Count the number of cells in the region and the ;argest region is the
% bridge bottom surface
num_region = max(largest_region_cell_ids);
region_count = histcounts(largest_region_cell_ids, num_region);
region_stat(:, 1) = unique(largest_region_cell_ids);
region_stat(:, 2) = region_count(region_count> 0);
[~, sort_ids] = sort(region_stat(:, 2), "descend");
bridge_bottom_surf_id = region_stat(sort_ids(1),1);
clear num_region region_count region_stat ids

%% Find the region connect to the reference surface
% component_link_list = [master id, slave id, type connection[-1, 0, 1], side connection[-1,1]; -1 coincide; 0 parallel; 1 intersect 
connect_segs_ids = setdiff(region_features.ids, bridge_bottom_surf_id); 
comp_link_list = connected_2D_components(Tree, cells_region, region_features,...
                        bridge_bottom_surf_id, connect_segs_ids, struct_threshold, threshold);
                    
un_connect_segs_ids = setdiff(region_features.ids,reshape(comp_link_list(:,[1,2]),[],1));

%% Filtering outlier points
mask = (comp_link_list(:,1) == bridge_bottom_surf_id)&(comp_link_list(:,3) == 1);
filter_regions_ids = comp_link_list(mask,[1,2]);
for i = 1:size(filter_regions_ids,1)
    cells_region = plane_outlier_removal(Tree, cells_region, filter_regions_ids(i,:), threshold);
end

%% Assign comp_link_list to the bridge component
% Update bridge tree
% For bottom surface of the superstructure
bridge_bottom_name = 'Bottom Beam';
[supstruct, cells_region] = bridge_tree_components(supstruct, ' ', bridge_bottom_surf_id, bridge_bottom_name, cells_region);


%% Assign comp_link_list to the bridge component
master_comp_ids = unique(comp_link_list(:,1));
master_name_list = cell(numel(master_comp_ids),1);
mask = ismember(master_comp_ids,bridge_bottom_surf_id);
master_name_list{mask} = bridge_bottom_name;
flag = true(size(comp_link_list,1),1);
while any(flag)
    % Retrieve the connection
    connect_id = find(flag, 1, 'first');
    connect_list = comp_link_list(connect_id,:);
    
    % Master segment
    master_comp_id = connect_list(1);
    mask = ismember(master_comp_ids, master_comp_id);
    master_comp_name = master_name_list{mask};
    
    % Slave segment
    % Retrieve all segments connected to the master on the same side
    mask = ismember(comp_link_list(:,1), connect_list(1)) & ismember(comp_link_list(:,4), connect_list(4));
    slave_comp_ids = comp_link_list(mask,2);
    
    % Get name for the slave
    slave_comp_name = master_comp_name;

    % Add a connection side   
    if connect_list(4) == 1
        slave_comp_name = strcat(slave_comp_name, " ", num2str(1));
    else
        slave_comp_name = strcat(slave_comp_name, " ", num2str(2));
    end
    
    % Assign the data structure
    [supstruct, cells_region] = bridge_tree_components(supstruct, master_comp_name, slave_comp_ids, slave_comp_name, cells_region);
    
    % Mask segment are assigned
    mask = ismember(comp_link_list(:,1), master_comp_id)&ismember(comp_link_list(:,2), slave_comp_ids);
    flag(mask) = false;
    
    % Assign the slave name list
    mask = ismember(master_comp_ids, slave_comp_ids);
    if any(mask)
        master_name_list{mask} = slave_comp_name;
    end
end



