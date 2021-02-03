function [supstruct, cells_region, un_connect_segs_ids, road_surface] = bridge_top_superstructure(Tree, cells_region, supstruct, struct_threshold, threshold, time_record)
%% This function is to extract the components surface of the bridge structure
% Input:
%     Tree                      : Data structure stored 2D cell in xy associated with ptc_ids
%     supstruct                 : Data structure store results
%     cells_region              : Data structure stored points in each segment (plane)
%     struct_threshold          : Structur_threshold relating a minimum size of the structural components
%     threshold                 : Threashold relating to filtering and rocessing point cloud
%     debug                     : Check running time

% Output
%     supstruct                 : Data structure store results
%     cells_region              : Data structure stored points in each segment (plane)
%     un_connect_segs_ids       : [Nx1] list of non-segment ids
%     road_surface              : Data structure store road features

%% Demo
% Tree = OQTR;
% supstruct = SuperStructure;
% cells_region = Cells_Region;
% struct_threshold;
% threshold = THRESHOLD;
% time_record = true;

% Return results
% SuperStructure =supstruct;
% Cells_Region = cells_region;
% un_connect_segs_ids = un_connect_segs_ids;
% Road_Surface = road_surface;
%% Compute the feature of the region
region_features = cal_region_features(Tree, cells_region);
%% Find the reference component: top: road surface; 
% Road surface
[~, max_id] = max(region_features.bounding_box(:,1));
road_surf_id = region_features.ids(max_id, 1);
%     ref_seg_cell_ids = vertcat(cells_region(ref_seg_ids).cell.id);

% Update bridge tree
road_surf_name = 'Road Surface';
[supstruct, cells_region] = bridge_tree_components(supstruct, ' ', road_surf_id, road_surf_name, cells_region);

%% Find the region connect to the reference surface
mask = (region_features.center(:,3) >= region_features.center(road_surf_id,3))&...
       (region_features.bounding_box(:,2) >= 0.5*struct_threshold.bridge_min_length);
connect_segs_ids = region_features.ids(mask,1);
connect_segs_ids = setdiff(connect_segs_ids, road_surf_id); 
un_connect_segs_ids = region_features.ids(~mask,1);
comp_link_list = connected_2D_components(Tree, cells_region, region_features,...
                        road_surf_id, connect_segs_ids, struct_threshold, threshold);

% Check if child on the same side and connect together
% flag = true(size(comp_link_list,1),1);
% for i = 1:size(comp_link_list,1)
%     connect_comp_ids = comp_link_list(i,:);
%     mask = ismember(comp_link_list(:,[1,2]),connect_comp_ids(1:2), 'rows');
%     remain_comp_link_list = comp_link_list(~mask,:);
%     mask = ismember(remain_comp_link_list(:,2), reshape(connect_comp_ids(1:2),[],1));
%     parent_id = (remain_comp_link_list(mask, 1));
%     side_direction = (remain_comp_link_list(mask, 3));
%     if (numel(parent_id) == 2)&&(numel(unique(parent_id)) == 1)&&(numel(side_direction) == 2)&&(numel(unique(side_direction) == 1))
%         flag(i) = false;
%     end
% end
% comp_link_list = comp_link_list(flag,:);
% clear i mask connect_component_ids remain_component_link_list parent_id side_direction flag

%% Assign comp_link_list to the bridge component
master_comp_ids = unique(comp_link_list(:,1));
master_name_list = cell(numel(master_comp_ids),1);
mask = ismember(master_comp_ids, road_surf_id);
master_name_list{mask} = road_surf_name;
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
    % Get name for the slave
    if master_comp_id == road_surf_id
        slave_comp_name = 'Footpath';
    else
        slave_comp_name = master_comp_name;
    end
    
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

%% Compute features of the road surface
road_surface = struct('ptc_ids', [], 'cent', [], 'normal',[], 'tangent',[], 'mbb',[]);
mask = cellfun(@(s) contains('Road Surface', s, 'IgnoreCase', true), supstruct(1).Description);
road_surf_ptc_ids = vertcat(supstruct.Component(mask).cell.ptc_ids);
road_surf_ptc_xyz = Tree.pts(road_surf_ptc_ids,:);
road_surf_eigenspace = eigenspace(road_surf_ptc_xyz, 0);

% Estimate the bounding box
rot_euler_angle = vrrotvec(road_surf_eigenspace(1,:),threshold.nz);
rot_matrix = vrrotvec2mat(rot_euler_angle);
rot_road_surf_ptc_xyz = (rot_matrix*road_surf_ptc_xyz')';
rot_road_mbb = min2DBoundingBox(rot_road_surf_ptc_xyz(:,1:2)');

% Assigned to the structure    
road_surface.ptc_ids = road_surf_ptc_ids;
road_surface.cent = mean(road_surf_ptc_xyz,1);
road_surface.normal = road_surf_eigenspace(1,:);
road_surface.tangent = road_surf_eigenspace(3,:);

% Update thrid direction for mmb
rot_road_mbb.vertices(:,3) = road_surface.cent(3);
rot_road_mbb.long_edge_vector(:,3) = 0;
rot_road_mbb.short_edge_vector(:,3) = 0;
rot_road_mbb.long_edge_vector = rot_road_mbb.long_edge_vector/norm(rot_road_mbb.long_edge_vector);
rot_road_mbb.short_edge_vector = rot_road_mbb.short_edge_vector/norm(rot_road_mbb.short_edge_vector);
road_surface.mbb = rot_road_mbb;


   