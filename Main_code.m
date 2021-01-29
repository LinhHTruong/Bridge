%% Read file
clc; close all; clear;
curdir = '\Concrete Bridge\';
filename = 'Concrete_Bridge_10mm.txt';

delimiter = {' ' '\b' '\t' ',' ';' }; % all possible delimiter 
fid = fopen(strcat(curdir,filename));
line = fgetl(fid);
firstLine = textscan(line, '%f', 'Delimiter', delimiter, 'MultipleDelimsAsOne', 1);
nocolumn = length(firstLine{1,1});
formatSpec = [repmat('%f ', 1, nocolumn) '%*[^\n]'];
frewind(fid); % jump back to beginning of file
% Read all data
allData = textscan(fid,formatSpec,'Delimiter', delimiter,'MultipleDelimsAsOne', 1);
fclose(fid);
%
PTC.xyz = cell2mat(allData);
clear allData
% PTC.xyz = PTC.xyz(:,1:3);
clear fid line firstLine nocolumn delimiter formatSpec %allData
%% Estabish threshols
sampling_step = 0.005;      % m, Real sampling step of the data
cell_size = 0.5;            % m, Input cell size and will be adjust with real cell size
voxel_size = 0.05;          % m, Input voxel size
max_angle = 5.0;            % degree, Maximum deviation angle between two surfaces (local or global)
max_hor_plane_angle = 60;   % degree, The maximum angle between the horizontal plane and vertical
max_dev_dist = 0.01;        % m, A maximum deviation distance between the local surfaces
dev_plane = 0.01;           % Deviation 
residual = 0.005;           % m, Residual of the (local) surface 
distance = 0.01;
min_num_pts = 10;
min_num_cell = 3;
max_slope = 0.1;            % Maximum a slope

THRESHOLD = Sematic_thresholds('sampling_step',sampling_step,'voxel_size',voxel_size,'cell_size',cell_size, 'max_angle', max_angle, 'max_slope',max_slope,...
    'max_dev_dist', max_dev_dist,'dev_plane',dev_plane, 'max_hor_plane_angle',max_hor_plane_angle, 'min_num_pts', min_num_pts, 'residual', residual,'distance', distance,'min_num_cell',min_num_cell);
clear cell_size voxel_size max_angle max_slope max_height dev_height min_num_pts residual distance min_num_cell max_dev_angle_hor_plane max_dev_height   
% Generate cell grids
tic
OQTR = OctQuadtree(PTC.xyz,'max_size', THRESHOLD.cell_size);
fprintf('Running time of generating quadtree: %.2f seconds \n', toc);
% clear PTC
% Create data structure stored point clouds of the bridge components
Bridge = struct('Component',[], 'Link_List',[]);
%% Extract the cells on a leaf node containing the horizontal plane
leaf_cell_ids = Node_Leaf(OQTR);
leaf_cell_ids = leaf_cell_ids(OQTR.cell_props(leaf_cell_ids) == 1);
smallest_cell_size = (OQTR.cell_bounds(leaf_cell_ids,[4,5]) - OQTR.cell_bounds(leaf_cell_ids,[1,2]));
THRESHOLD.cell_size = max(smallest_cell_size(:));
[Plane_Cells, nonPlane_Cells] = Cell_Planes_Extraction(OQTR, leaf_cell_ids, THRESHOLD);
clear smallest_cell_size
%% Region growing
[Region_Info, Plane_Cells] = bridge_cell_segmentation(OQTR, Plane_Cells, THRESHOLD, 'top');
%% Back and Forward filtering: For the cells are shared by multiple segments/regions
Region = backward_forward_filtering(OQTR, Plane_Cells, Region_Info, THRESHOLD);
%
%% Forward region growing: for the cells on the boundary of the region to
% searching points out off the peak (the surface) of the cells out off the region

[Plane_Cells, Region] = forward_cell_region_growing(OQTR, Plane_Cells, Region, THRESHOLD);

%% Back and Forward filtering: For the cells are shared by multiple segments/regions
[Plane_Cells, Region] = boundary_filtering(OQTR, Plane_Cells, Region, THRESHOLD);

%% Forward searching the points on non-planar cells
[nonPlane_Cells, Region] = non_plane_cell_filtering(OQTR, nonPlane_Cells, Region, THRESHOLD);
%% Merging regions
Region = merge_segments(OQTR, Region, THRESHOLD);

%% Compute features of each region
tic
Region_Features = cal_region_features(Region, OQTR);
fprintf('Running time of computing features of the regions: %.2f seconds \n', toc);
%% Built the data structure to identify the bridge component

bridge_threshold.footpath_min_width = 0.8;
bridge_threshold.footpath_length_ratio = 0.25;
bridge_threshold.parapet_min_height = 0.6;
bridge_threshold.curb_max_height = 0.3; % Maximum height of the road curb, which mean the different height between the road surface and footpath is less than this threshold
bridge_threshold.column_min_height = 0.5;
bridge_threshold.adjoin_components_min_edge_length = 2.0; % The minimum length of edges sharing by two adjoin or parallel surfaces
bridge_threshold.surface_error = 0.005;
bridge_threshold.parapet_rail_cluster = 0.2;
bridge_threshold.min_no_cells = 20;
bridge_threshold.bandwidth = 0.1;
bridge_threshold.search_scale = 2; % Window size to search neighbour cells
% 
% bridge_threshold.surface_normal_deviation = sin(deg2rad(5)); % The deviation of normals between two surfaces
% bridge_threshold.cell_size = THRESHOLD.cell_size;

bridge_tree_level = 1;
[Bridge, un_components_regions_ids] = bridge_desk_components(OQTR, Region, Region_Features, Bridge, bridge_tree_level, THRESHOLD, bridge_threshold);

% Update the status for peaks in the cells within the Plane_Cells
tic
Plane_Cells = region_status_update(OQTR, Plane_Cells, Region, Bridge, bridge_tree_level, un_components_regions_ids, THRESHOLD);
fprintf('Running time of update status of cells: %.2f seconds \n', toc);
%% Building up data structure for second level
tic
bridge_tree_level = 2;
[Bridge, un_components_regions_ids] = bridge_bottom_desk_components(OQTR, Plane_Cells, Region, Region_Features, ...
                                                Bridge, bridge_tree_level, THRESHOLD, bridge_threshold);

Plane_Cells = region_status_update(OQTR, Plane_Cells, Region, Bridge, bridge_tree_level, un_components_regions_ids, THRESHOLD);

fprintf('Running time of assigning components of the bottom surface of the bridge: %.2f seconds \n', toc);


%% For ground
% [Region_Info, Plane_Cells] = bridge_cell_segmentation(OQTR, Plane_Cells, THRESHOLD, 'bottom');
% %%
% Region = backward_forward_filtering(OQTR, Plane_Cells, Region_Info, THRESHOLD);
% %
% %% Forward region growing: for the cells on the boundary of the region to
% % searching points out off the peak (the surface) of the cells out off the region
% [Plane_Cells, Region] = forward_cell_region_growing(OQTR, Plane_Cells, Region, THRESHOLD);
% 
% %% Back and Forward filtering: For the cells are shared by multiple segments/regions
% [Plane_Cells, Region] = boundary_filtering(OQTR, Plane_Cells, Region, THRESHOLD);
% %% Forward searching the points on non-planar cells
% [nonPlane_Cells, Region] = non_plane_cell_filtering(OQTR, nonPlane_Cells, Region, THRESHOLD);
% %% Merging regions
% Region = merge_segments(OQTR, Region, THRESHOLD);
% 
% %% Update bridge tree for ground
% bridge_tree_level = 3;
% for i=1:length(Region)
%     [Bridge,~] = bridge_tree_components(Bridge, bridge_tree_level, 0, i, Region);
% end
%% Update a new cell points for the substructure
tic
bridge_substructure = struct('cell',[]);
for i = 1:numel(leaf_cell_ids)
    leaf_cell_id = leaf_cell_ids(i);
    bridge_substructure.cell(i).ids = leaf_cell_id;
    bridge_substructure.cell(i).ptc_ids = OQTR.cell_pts(leaf_cell_id).id;
end

% Filtering
for i = 1:length(Bridge)
    num_components = length(Bridge(i).Component);
    for j = 1:num_components
        num_cells = length(Bridge(i).Component(j).cell);
        for k = 1:num_cells
            % Retrieve the cell and points in the bridge components
            cell_id = Bridge(i).Component(j).cell(k).ids(1);
            cell_ptc_ids = Bridge(i).Component(j).cell(k).ptc_ids;
            
            % Update points
            mask = ismember(leaf_cell_ids, cell_id);
            substructure_cell_ptc_ids = bridge_substructure.cell(mask).ptc_ids;
            bridge_substructure.cell(mask).ptc_ids = setdiff(substructure_cell_ptc_ids,cell_ptc_ids);        
        end
    end
end

% Remove empty cells
temp = struct('cell',[]);
count = 1;
for i = 1:length(bridge_substructure.cell)
    cell_ptc_ids = bridge_substructure.cell(i).ptc_ids;
    if numel(cell_ptc_ids) >= THRESHOLD.min_num_pts
        temp.cell(count).ids = bridge_substructure.cell(i).ids;
        temp.cell(count).ptc_ids = cell_ptc_ids;
        count = count + 1;
    end
end
bridge_substructure = temp;
bridge_substructure_cell_ids = vertcat(bridge_substructure.cell.ids);
clear i j k temp cell_id cell_ptc_ids substructure_cell_ptc_ids
toc

%% Substracture - Column extraction
% Calculate height of the cell
tic
cell_height_info = inf(numel(bridge_substructure_cell_ids), 3); %[cell_id, largest_region_length]
cell_height_info(:,1) = bridge_substructure_cell_ids;
for i = 1:numel(bridge_substructure_cell_ids)
    % Compute height of the cell
    substructure_cell_ptc_ids = bridge_substructure.cell(i).ptc_ids;
    substructure_cell_ptc_z = OQTR.pts(substructure_cell_ptc_ids,3);
    substructure_cell_ptc_z = sort(substructure_cell_ptc_z);
    cell_height_info(i,2) = max(diff(substructure_cell_ptc_z));
    cell_height_info(i,3) = abs(substructure_cell_ptc_z(1) - substructure_cell_ptc_z(end));
end
mask = isinf(cell_height_info(:,2));
cell_height_info = cell_height_info(~mask,:);
fprintf('Running time of calculating ceil height: %.2f seconds \n', toc);
clear i mask substructure_cell_ptc_z substructure_cell_ptc_ids
%% Extract the supstructurce cells
% Filtering based on a cell height
mask = (bridge_threshold.column_min_height <= cell_height_info(:,3));
substructure_cell_info = cell_height_info(mask,:);
non_substructure_cell_info = cell_height_info(~mask,:);
substructure_region_cell_info = cells_connectivity(OQTR, substructure_cell_info(:,1));

%% Filtering 
% Extract the cell containing points of the bottom bridge deck
components_link_list = Bridge(2).Link_List;
bridge_bottom_deck_component_ids = components_link_list(components_link_list(:,1) == 0,2);
bridge_bottom_deck_cell_ids = (arrayfun(@(x) vertcat(Bridge(2).Component(x).cell.ids),bridge_bottom_deck_component_ids,'UniformOutput',false));
bridge_bottom_deck_cell_ids = cell2mat(bridge_bottom_deck_cell_ids);
clear components_link_list bridge_bottom_deck_component_ids

% Check cells in the region      
region_ids = unique(substructure_region_cell_info(:,2));
region_score = zeros(numel(region_ids),2);
region_score(:,1) = region_ids;
for i = 1:numel(region_ids)
    % Retrieve the points within the region
    mask = substructure_region_cell_info(:,2) == region_ids(i);
    region_cell_ids = substructure_region_cell_info(mask,1);
    mask = ismember(bridge_substructure_cell_ids, region_cell_ids);
    region_cell_ptc_ids = vertcat(bridge_substructure.cell(mask).ptc_ids);
    region_cell_ptc_xyz = OQTR.pts(region_cell_ptc_ids,1:3);
    min_boudingbox = min2DBoundingBox(region_cell_ptc_xyz(:,1:2)');
    
    % Extract cell on a boundary of the region: region having boundary +
    % long edge >3.0 me = 1 lane or 1/2 bridge width
    [region_bound_cell_ids, ~] = boundary_cells(OQTR, region_cell_ids);

    if (~isempty(region_bound_cell_ids))&&(min_boudingbox.long_edge >= 3.0)
        
        % Preallocate region_bound_cells
        region_bound_cell_info = inf(numel(region_bound_cell_ids),3);% [cell id, cell length, cell gap]
        region_bound_cell_info(:,1) = region_bound_cell_ids;
        for j = 1:numel(region_bound_cell_ids)
            % Check the highest point in the cell to the bottom surface of the slab
            mask = ismember(bridge_substructure_cell_ids, region_bound_cell_ids(j));
            cell_ptc_ids = bridge_substructure.cell(mask).ptc_ids;
            cell_ptc_z = OQTR.pts(cell_ptc_ids,3);

            % Find neighbour cells 
            adjacent_cell_ids = Query_Neighbour_Cells(OQTR, region_bound_cell_ids(j), bridge_bottom_deck_cell_ids(:,1));
            mask = ismember(bridge_bottom_deck_cell_ids, adjacent_cell_ids);
            adjacent_cell_ids = bridge_bottom_deck_cell_ids(mask,:);
            
            % Compute features of the cell: cell_height and gap between the cell and the bottom surface of the bridge deck
            if ~isempty(adjacent_cell_ids)
                mask = ismember(Plane_Cells.cell_ids(:,[1,2]),adjacent_cell_ids, 'rows');
                adjacent_cell_peak_features = vertcat(Plane_Cells.peak_info(mask).peaks_features);
                adjacent_cell_center_z = mean(adjacent_cell_peak_features(:,3));
                % Different between the highest points and the bottom surface of the bridge slab
                region_bound_cell_info(j,2:3) = [max(cell_ptc_z) - min(cell_ptc_z), abs(max(cell_ptc_z) - adjacent_cell_center_z)];
            end
        end
        
        % Build up score matrix
        mask = region_bound_cell_info(:,3) <= 0.5;
        region_bound_cell_info = region_bound_cell_info(mask,:);
        if size(region_bound_cell_info,1)*THRESHOLD.cell_size >= 3.0 %Update by a half width of the bridge 
            region_score(i,2) = 1;
        end
        clear j mask region_bound_cell_info cell_ptc_ids cell_ptc_z adjacent_cell_ids adjacent_cell_peak_features adjacent_cell_center_z
    end
end

% Update supstructure
substructure_component_region_ids = region_score(logical(region_score(:,2)),1);
mask = ismember(substructure_region_cell_info(:,2), substructure_component_region_ids);
substructure_region_cell_info = substructure_region_cell_info(mask,:);
clear i mask region_ids region_cell_ids region_cell_ptc_ids region_cell_ptc_xyz min_boudingbox region_bound_cell_ids

%% Classify the subsubstructure: Piers and Abutments
substructure_ids = unique(substructure_region_cell_info(:,2));
substructures_features = zeros(numel(substructure_ids),7);%[id, cent_x, cent_y, cent_z, long_edge, short_edge, height]
substructures_features(:,1) = substructure_ids;
for i = 1:numel(substructure_ids)
    mask = substructure_region_cell_info(:,2) == substructure_ids(i);
    substructure_cell_ids = substructure_region_cell_info(mask,1);
    mask = ismember(bridge_substructure_cell_ids, substructure_cell_ids);
    substructure_cell_ptc_ids = vertcat(bridge_substructure.cell(mask).ptc_ids);
    substructure_cell_ptc_xyz = OQTR.pts(substructure_cell_ptc_ids,1:3);

    substructure_center = mean(substructure_cell_ptc_xyz, 1);
    min_boudingbox = min2DBoundingBox(substructure_cell_ptc_xyz(:,1:2)');
    substructures_features(i,2:end) = [substructure_center, min_boudingbox.long_edge,min_boudingbox.short_edge, max(substructure_cell_ptc_xyz(:,3)) - min(substructure_cell_ptc_xyz(:,3))];
    clear substructure_cell_ptc_ids substructure_cell_ptc_xyz substructure_center min_boudingbox
end
fprintf('Running time of computing features of the regions: %.2f seconds \n', toc);
clear i mask substructure_ids substructure_cell_ptc_ids  substructure_cell_ptc_xyz substructure_center min_boudingbox

% Determine id of the substructure as abutments (end points) orther as piers
[~, ~, ~, start_ptc_id, end_ptc_id] = findEndPts(substructures_features(:,2:4));
pier_ids = setdiff([1:size(substructures_features,1)]',[start_ptc_id,end_ptc_id]);
abutment_ids = [start_ptc_id, end_ptc_id];
clear i mask substructure_ids substructures_features start_ptc_id end_ptc_id 
%% Abutment segmentation
bridge_tree_level = length(Bridge);
for i = 1:numel(abutment_ids)    
    
    % Retrieve the points of the abutment
    mask = substructure_region_cell_info(:,2) == abutment_ids(i);
    substructure_cell_ids = substructure_region_cell_info(mask,1);
    mask = ismember(bridge_substructure_cell_ids, substructure_cell_ids);
    substructure_cell_ptc_ids = vertcat(bridge_substructure.cell(mask).ptc_ids);
%     substructure_cell_ptc_xyz = OQTR.pts(substructure_cell_ptc_ids,1:3);
    
    % call clustering function here :)
    Bridge = abutment_segmentation(OQTR, substructure_cell_ptc_ids, Bridge, bridge_tree_level + i, THRESHOLD);

end

%% Pier segmentation
bridge_tree_level = length(Bridge);
for i = 1:numel(pier_ids)    
    
    % Retrieve the points of the abutment
    mask = substructure_region_cell_info(:,2) == pier_ids(i);
    substructure_cell_ids = substructure_region_cell_info(mask,1);
    mask = ismember(bridge_substructure_cell_ids, substructure_cell_ids);
    substructure_cell_ptc_ids = vertcat(bridge_substructure.cell(mask).ptc_ids);
%     substructure_cell_ptc_xyz = OQTR.pts(substructure_cell_ptc_ids,1:3);
    
    % call clustering function here :)
%     Bridge = abutment_segmentation(OQTR, substructure_cell_ptc_ids, Bridge, bridge_tree_level + i, THRESHOLD);

end

























































%% Column
cur_horizon_level = length(Bridge_Component);

surface_ptc_local_ids = vertcat(Pier.ptc_ids);   
surface_ptc_global_ids = pier_cell_ptc_ids(surface_ptc_local_ids);
surface_ptc_xyz = OQTR.pts(surface_ptc_global_ids,1:3);
    

Bridge_Component(cur_horizon_level+1).Description = 'column';

Bridge_Component(cur_horizon_level+1).cell.ptc_ids = surface_ptc_global_ids;
Bridge_Component(cur_horizon_level+1).cell.ptc_xyz = surface_ptc_xyz;



%% Combinate all points
classification_ptc = zeros(size(OQTR.pts,1),1);
cur_horizon_level = length(Bridge_Component);
for i = 1:cur_horizon_level
   level_ptc_ids =  vertcat(Bridge_Component(i).cell.ptc_ids);
   classification_ptc(level_ptc_ids) = i;
end
ground_ptc_ids = vertcat(Ground.cell.ptc_ids);
classification_ptc(ground_ptc_ids) = max(classification_ptc) + 1;
%% Write data
outputFile = strcat(our_dir, 'Germany_Concrete_Ptc_Classification.txt');
write_txt(outputFile, [OQTR.pts,classification_ptc])








%% Write the cell region
all_ptc_xyz_class = [];
for i = 1:max(substructure_region_cell_info(:,2))
    mask = substructure_region_cell_info(:,2) == i;
    region_cell_ids = substructure_region_cell_info(mask,1);
    region_cell_ptc_ids = vertcat(OQTR.cell_pts(region_cell_ids).id);
    region_cell_ptc_xyz = OQTR.pts(region_cell_ptc_ids,1:3);
    region_cell_ptc_xyz(:,4) = i;
    all_ptc_xyz_class = [all_ptc_xyz_class;region_cell_ptc_xyz]; 
end







%% recover un-classification cell 
non_column_cell_center = (OQTR.cell_bounds(non_substructure_cell_info(:,1),1:2)+...
                          OQTR.cell_bounds(non_substructure_cell_info(:,1),4:5))/2.0;
             
for i = 1:numel(bridge_support_component_region_ids)
    % Retrieve cells in the region
    mask = substructure_region_cell_info(:,2) == bridge_support_component_region_ids(i);
    region_cell_ids = substructure_region_cell_info(mask,1);
    region_cell_ptc_ids = vertcat(OQTR.cell_pts(region_cell_ids).id);
    region_cell_ptc_xyz = OQTR.pts(region_cell_ptc_ids,1:3);
    min_boudingbox = min2DBoundingBox(region_cell_ptc_xyz(:,1:2)');
    min_boudingbox_polygon = [min_boudingbox.vertices;min_boudingbox.vertices(1,:)];
    
    % Check cells in 
    in = inpolygon(non_column_cell_center(:,1),non_column_cell_center(:,2),...
        min_boudingbox_polygon(:,1),min_boudingbox_polygon(:,2));
    
    add_column_cell_ids = non_substructure_cell_info(in,1);
    if ~isempty(add_column_cell_ids)
        add_column_cell_ids(:,2) = bridge_support_component_region_ids(i);
        substructure_region_cell_info = union(substructure_region_cell_info, add_column_cell_ids, 'rows');
        non_column_cell_center(in,:) = [];
        non_substructure_cell_info(in,:) = [];
    end
end

%% Clean up the points above the level
% OQTR = cleanning_cells(Bridge, bridge_tree_level, OQTR, THRESHOLD, bridge_threshold.surface_error, 'lowest');
% 
% 
% 
% %%
    
% elseif horizontal_level == 2 % Bottom surface of the bridge girder
%     % Determine the bottom surface of the slab/girder: large regions,
%     % coplanar
%     [~, max_id] = max(Region_Features(:,11));
%     bridge_bottom_desk_region_ids = Region_Features(max_id, 1);
%     bridge_bottom_desk_region_cent = Region_Features(max_id, 2:4);
%     bridge_bottom_desk_region_normal = Region_Features(max_id, 5:7);
%     bridge_bottom_desk_region_tangent = Region_Features(max_id,8:10);
%     bridge_bottom_desk_region_bound = Region_Features(max_id, 12:14);
%     
%     % Find co-planar region: cosin of normals and orthogonal distance;
%     % similar in geometry: width no large different
%     cosine_normals = cosine_vectors(bridge_bottom_desk_region_normal, Region_Features(:, 5:7));
%     dist_region_2_region = distPts2Surf(Region_Features(:,2:4), bridge_bottom_desk_region_cent, bridge_bottom_desk_region_normal);
%     
%     mask = ((THRESHOLD.max_angle <= abs(cosine_normals))&...
%             (dist_region_2_region <= 20*THRESHOLD.dev_plane)&...
%             (0.75*bridge_bottom_desk_region_bound(3)<=Region_Features(:, 14)));
%     bridge_structure_region_ids = Region_Features(mask,1); % Inlude the largest region
%     bridge_structure_region_cent = Region_Features(mask, 2:4);
%     bridge_structure_region_bound = Region_Features(mask, 12:14);
%     % Assign the bridge components: bottom surfaces of the bridge desk
%     cur_horizon_level = length(Bridge_Component);
%     for i = 1:numel(bridge_structure_region_ids)
%         Bridge_Component(cur_horizon_level+i).Description = strcat('bridge slab',{' '},num2str(i));
%         Bridge_Component(cur_horizon_level+i).cell = Region(bridge_structure_region_ids(i)).cell;
%         Bridge_Component(cur_horizon_level+i).features = [Region_Features(bridge_structure_region_ids(i),2:4),...
%                                                           Region_Features(bridge_structure_region_ids(i),8:10),...
%                                                           Region_Features(bridge_structure_region_ids(i),12:14)];
%     end
%     clear cosine_normals dist_region_2_region mask i
%     % Update region features
%     mask = ismember(Region_Features(:,1),bridge_structure_region_ids);
%     Region_Features(mask,:) = [];
%     
%     % Update bottom surface
%     bridge_supstructure_cell_ptc_xyz = (arrayfun(@(x) vertcat(Region(x).cell.ptc_xyz),bridge_structure_region_ids,'UniformOutput',false));
%     bridge_supstructure_cell_ptc_xyz = cell2mat(bridge_supstructure_cell_ptc_xyz);
%     min_boudingbox = min2DBoundingBox(bridge_supstructure_cell_ptc_xyz(:,1:2)');
%     bridge_bottom_slab_center = mean(bridge_supstructure_cell_ptc_xyz);
%     region_eigenspace = eigenspace(bridge_supstructure_cell_ptc_xyz, 0);
%     bridge_bottom_slab_tangent = region_eigenspace(3,:);
%     bridge_bottom_slab_width = min_boudingbox.short_edge;
%     bridge_bottom_slab_length = min_boudingbox.long_edge;
%     clear bridge_supstructure_cell_ptc_xyz min_boudingbox region_eigenspace
%     
%     % Extract the large regions
%     mask = (bridge_bottom_slab_center(3) < Region_Features(:, 4))&...
%            (0.5*bridge_bottom_slab_length <= Region_Features(:,13));
%     side_bottom_slab_regions_ids = Region_Features(mask,1);
%     if ~isempty(side_bottom_slab_regions_ids)
%         side_bottom_slab_regions_center = Region_Features(mask,2:4);
%         side_bottom_slab_regions_tangent = Region_Features(mask,8:10);
% 
% 
%         % Check footpath region on the same side
%         side_bottom_slab_cent_sign = ones(numel(side_bottom_slab_regions_ids),1);
%         side_bottom_slab_proj_cent = proj_point_on_3Dline(side_bottom_slab_regions_center, [bridge_bottom_slab_center, bridge_bottom_slab_tangent]);
%         side_bottom_slab_cent_vect = side_bottom_slab_regions_center - side_bottom_slab_proj_cent;
%         side_bottom_slab_cent_sign(2:end) = sign(sum(bsxfun(@times,side_bottom_slab_cent_vect(1,:),side_bottom_slab_cent_vect(2:end,:)),2));
% 
%         % Assign the bridge components: footpath
%         cur_horizon_level = length(Bridge_Component);
%         side_bottom_slab_cent_group = unique(side_bottom_slab_cent_sign);
%         for i = 1:numel(side_bottom_slab_cent_group)%numel(footpath_regions_ids)
%             group_id = side_bottom_slab_cent_sign(i);
%             mask = side_bottom_slab_cent_sign == group_id;
%             side_bottom_slab_regions_id = side_bottom_slab_regions_ids(mask);
% 
%             % Update bridge component
%             count = 1;
%             Bridge_Component(cur_horizon_level+i).Description = strcat('bridge side slab',{' '},num2str(i));
%             for j = 1:numel(side_bottom_slab_regions_id)
%                region_length = length(Region(side_bottom_slab_regions_id(j)).cell);
%                for k = 1:region_length
%                    Bridge_Component(cur_horizon_level+i).cell(count).id = Region(side_bottom_slab_regions_id(j)).cell(k).id;
%                    Bridge_Component(cur_horizon_level+i).cell(count).ptc_ids = Region(side_bottom_slab_regions_id(j)).cell(k).ptc_ids;
%                    Bridge_Component(cur_horizon_level+i).cell(count).ptc_xyz = Region(side_bottom_slab_regions_id(j)).cell(k).ptc_xyz;
%                    Bridge_Component(cur_horizon_level+i).cell(count).surface_centroid = Region(side_bottom_slab_regions_id(j)).cell(k).surface_centroid;
%                    Bridge_Component(cur_horizon_level+i).cell(count).surface_normls = Region(side_bottom_slab_regions_id(j)).cell(k).surface_normls;
%                    Bridge_Component(cur_horizon_level+i).cell(count).status = Region(side_bottom_slab_regions_id(j)).cell(k).status;
% 
%                    count = count + 1;
%                end 
%             end
%             % Update feature
%             if numel(side_bottom_slab_regions_id) >= 2
% 
%                 region_ptc_xyz = (arrayfun(@(x) vertcat(Region(x).cell.ptc_xyz),side_bottom_slab_regions_id,'UniformOutput',false));
%                 region_ptc_xyz = cell2mat(region_ptc_xyz);
%                 region_cent = mean(region_ptc_xyz, 1);
%                 region_eigenspace = eigenspace(region_ptc_xyz, 0);
%                 region_tangent = region_eigenspace(3,:);
%                 min_boudingbox = min2DBoundingBox(region_ptc_xyz(:,1:2)');
%                 Bridge_Component(cur_horizon_level+i).features = [region_cent,region_tangent,min_boudingbox.long_edge,min_boudingbox.short_edge];
%                 clear region_ptc_xyz region_ptc_xyz min_boudingbox
%             else
%                 Bridge_Component(cur_horizon_level+i).features = [Region_Features(side_bottom_slab_regions_id,2:4),...
%                                                                   Region_Features(side_bottom_slab_regions_id,8:10),...
%                                                                   Region_Features(side_bottom_slab_regions_id,12:14)];
%             end
%             clear group_id mask footpath_regions_id
%         end
%         
%         % Update region feature
%         mask = ismember(Region_Features(:,1), side_bottom_slab_regions_ids);
%         Region_Features(mask,:) = [];
%     
%     end
%     
%     % Extract any cell
%     bridge_component_names = vertcat(Bridge_Component.Description);
%     bridge_desk_level_ids = find(contains(bridge_component_names, "slab")==0);
%     bottom_slab_level_ids = find(contains(bridge_component_names, "slab")==1);
%     
%     % retrieve cell is bridge desk+footpath +parapet
%     bridge_desk_cell_ids = (arrayfun(@(x) vertcat(Bridge_Component(x).cell.id),bridge_desk_level_ids,'UniformOutput',false));
%     bridge_desk_cell_ids = cell2mat(bridge_desk_cell_ids);
%     
%     bridge_desk_cell_surface_centroid = (arrayfun(@(x) vertcat(Bridge_Component(x).cell.surface_centroid),bridge_desk_level_ids,'UniformOutput',false));
%     bridge_desk_cell_surface_centroid = cell2mat(bridge_desk_cell_surface_centroid);
%     
%     % retrieve cell is bottom slab
%     bottom_slab_cell_ids = (arrayfun(@(x) vertcat(Bridge_Component(x).cell.id),bottom_slab_level_ids,'UniformOutput',false));
%     bottom_slab_cell_ids = cell2mat(bottom_slab_cell_ids);
%     bottom_slab_cell_surface_centroids = (arrayfun(@(x) vertcat(Bridge_Component(x).cell.surface_centroid),bottom_slab_level_ids,'UniformOutput',false));
%     bottom_slab_cell_surface_centroids = cell2mat(bottom_slab_cell_surface_centroids);
%     
% %     bottom_slab_cell_surface_normls = (arrayfun(@(x) vertcat(Bridge_Component(x).cell.surface_normls),bottom_slab_level_ids,'UniformOutput',false));
% %     bottom_slab_cell_surface_normls = cell2mat(bottom_slab_cell_surface_normls);
% %     cosine_val = cosine_vectors([0,0,1], bottom_slab_cell_surface_normls);
% %     mask = cosine_val < 0;
% %     bottom_slab_cell_surface_normls(mask,:) = -bottom_slab_cell_surface_normls(mask,:);
%     
%     % Extract the cell belonging to the top but not in the bottom slab
%     mask = ismember(bridge_desk_cell_ids,bottom_slab_cell_ids);% Cell in a desk + footpath+parapeter but not in bottom slaba 
%     bottom_slab_add_cell_ids = bridge_desk_cell_ids(~mask);
%     bottom_slab_add_cell_surface_centroids = bridge_desk_cell_surface_centroid(~mask,1:3);
%     
%     
%     % Clustering these cells -> divide location
%     add_slab = struct('cell',[]);
%     count = 1;
%     for i = 1:numel(bottom_slab_add_cell_ids)
%         % Determine the top surface threshold -> top surface: desk, parapet
%         bottom_slab_add_cell_id = bottom_slab_add_cell_ids(i);
%         if OQTR.cell_props(bottom_slab_add_cell_id)
%   
%             
%             bottom_slab_add_cell_surface_centroid = bottom_slab_add_cell_surface_centroids(i,:);
% 
%             % Find the closest cell -> bottom surface threshold: bottom slab
%             neighbour_cell_id = Closest_Cell(OQTR, bottom_slab_add_cell_id, bottom_slab_cell_ids);
%             mask = ismember(bottom_slab_cell_ids, neighbour_cell_id);
%             bottom_slab_cell_surface_centroid = mean(bottom_slab_cell_surface_centroids(mask,:),1);
%             
%             % Find the point within the cells
%             flag = true;
%             k = 1;
%             while (flag)&(k <= numel(bottom_slab_level_ids))
%                 temp_cell_ids = vertcat(Bridge_Component(bottom_slab_level_ids(k)).cell.id);
%                 temp_cell_id = find(temp_cell_ids == neighbour_cell_id);
%                 if ~isempty(temp_cell_id)
%                     bottom_slab_cell_ptc_min_z = min(vertcat(Bridge_Component(bottom_slab_level_ids(k)).cell(temp_cell_id).ptc_xyz(:,3)));
%                     flag = false;
%                 end
%                 k = k + 1;
%             end
%             
%             
% %             bottom_slab_cell_surface_norml = mean(bottom_slab_cell_surface_normls(mask,:),1);
% 
%             % Extract the points within the cell
%         
%             cell_ptc_ids = vertcat(OQTR.cell_pts(bottom_slab_add_cell_id).id);
%             cell_ptc_xyz = OQTR.pts(cell_ptc_ids,1:3);
% 
%             % Determine the points belonging to the slab: z from the bottom to
%             % the top
%             mask = (bottom_slab_cell_ptc_min_z <= cell_ptc_xyz(:,3))&(cell_ptc_xyz(:,3) <= bottom_slab_add_cell_surface_centroid(3));
%             slab_cell_ptc_ids = cell_ptc_ids(mask);
%             slab_cell_ptc_xyz = cell_ptc_xyz(mask,:);
% 
%             % Add new points for slab
%             if numel(slab_cell_ptc_ids) >= 20
%                 add_slab.cell(count).id = bottom_slab_add_cell_id;
%                 add_slab.cell(count).ptc_ids = slab_cell_ptc_ids;
%                 add_slab.cell(count).ptc_xyz = slab_cell_ptc_xyz;
%                 add_slab.cell(count).surface_centroid = mean(slab_cell_ptc_xyz,1);
%                 count = count + 1;
%             end
% 
%             % Update points for OQTR
%             mask = (cell_ptc_xyz(:,3) <= bottom_slab_cell_surface_centroid(3) - 0.01);
%             cell_sub_ptc_ids = cell_ptc_ids(mask);
%             OQTR.cell_pts(bottom_slab_add_cell_id).id = cell_sub_ptc_ids;
%             if numel(cell_sub_ptc_ids) > 20
%                 OQTR.cell_props(bottom_slab_add_cell_id) = 0;
%             end
%         end
%     end
%     
%     % Update bridge component: Slab Side
%     add_slab_ptc_xyz = vertcat(add_slab.cell.ptc_xyz);
%     min_boudingbox = min2DBoundingBox(add_slab_ptc_xyz(:,1:2)');
%   
%     cur_horizon_level = length(Bridge_Component);
%     Bridge_Component(cur_horizon_level+1).Description = 'bridge add slab';
%     Bridge_Component(cur_horizon_level+1).cell = add_slab.cell;
%     Bridge_Component(cur_horizon_level+1).features = [mean(add_slab_ptc_xyz, 1),min_boudingbox.area,min_boudingbox.long_edge,min_boudingbox.short_edge];
% 
%     clear add_slab add_slab_ptc_xyz min_boudingbox
%     
%     % Remove points in OQTR based bottom slab cell
%     
% %     leaf_cell_ids = Node_Leaf(OQTR);
% %     leaf_cell_ids = leaf_cell_ids(OQTR.cell_props(leaf_cell_ids) == 1);
% %     
% %     for i = 1: numel(leaf_cell_ids)
% %         leaf_cell_id = leaf_cell_ids(i);
% %         leaf_cell_ptc_ids = OQTR.cell_pts(leaf_cell_id).id;
% %         leaf_cell_ptc_xyz = OQTR.pts(leaf_cell_ptc_ids,1:3);
% %         
% %         mask = (leaf_cell_ptc_xyz(:,3) <= bridge_bottom_slab_center(3) - 0.025); % 0.025 offset
% %         bridge_sub_part_ptc_ids = leaf_cell_ptc_ids(mask);
% %         OQTR.cell_pts(leaf_cell_id).id = bridge_sub_part_ptc_ids; % Remaining ptc_ids 
% %         if isempty(bridge_sub_part_ptc_ids)
% %             OQTR.cell_props(leaf_cell_id) = 0;
% %         end
% % 
% %     end
%     
%     for i = 1:numel(bottom_slab_level_ids) % Throught each component
%         bottom_slab_level_id = bottom_slab_level_ids(i);
%         bottom_slab_cell_ids = vertcat(Bridge_Component(bottom_slab_level_id).cell.id);
% 
%         for j = 1:numel(bottom_slab_cell_ids) % Go through each cells
%             
%             % For cell belong to the bottom slab
%             bottom_slab_cell_id = bottom_slab_cell_ids(j);
%             bottom_slab_cell_surface_centroid = Bridge_Component(bottom_slab_level_id).cell(j).surface_centroid;
%             bottom_slab_cell_ptc_ids = vertcat(Bridge_Component(bottom_slab_level_id).cell(j).ptc_ids);
%             
%             
%             % Extract all points in cells in OQTR
%             leaf_cell_ptc_ids = OQTR.cell_pts(bottom_slab_cell_id).id;
%             leaf_cell_ptc_xyz = OQTR.pts(leaf_cell_ptc_ids,1:3);
%             
%             mask = ismember(leaf_cell_ptc_ids,bottom_slab_cell_ptc_ids);
%             cell_remain_ptc_ids = leaf_cell_ptc_ids(~mask);
%             cell_remain_ptc_xyz = leaf_cell_ptc_xyz(~mask,1:3);
% 
% 
% 
%             % Extract points belonging to a lower part -> Update OQTR
%             mask = cell_remain_ptc_xyz(:,3) < bottom_slab_cell_surface_centroid(3) - 0.075;
%             bridge_sub_part_ptc_ids = cell_remain_ptc_ids(mask);
%             OQTR.cell_pts(bottom_slab_cell_id).id = bridge_sub_part_ptc_ids; % Remaining ptc_ids 
%             if isempty(bridge_sub_part_ptc_ids)
%                 OQTR.cell_props(bottom_slab_cell_id) = 0;
%             end
%         end
%         
%     end
%     
%     
% %     %
% %     % Extract the points of slab
% %     bridge_slab = struct('cell',[]); 
% %     count = 1;
% %     cur_horizon_level = length(Bridge_Component);
% %     % For the cells of the bottom surfaces
% %     % Note: bridge_component_cell_surface_centroid(3) - 0.01: for remove
% %     % the points within the bottom surface of the bridge desk; 0.01: for
% %     % noise and bearing
% %     for i=1:cur_horizon_level
% %         if ~isempty(cell2mat(strfind(Bridge_Component(i).Description, 'bridge desk bottom')))
% %             bridge_component_cell_ids = vertcat(Bridge_Component(i).cell.id);
% % %             bridge_component_cell_surface_centroid = vertcat(Bridge_Component(i).cell.surface_centroid);
% %             
% %             for j = 1:numel(bridge_component_cell_ids)
% %                 % Retrieve points in bridge commponent cell
% %                 bridge_component_cell_ptc_ids = Bridge_Component(i).cell(j).ptc_ids;
% %                 bridge_component_cell_surface_centroid = vertcat(Bridge_Component(i).cell(j).surface_centroid);
% %             
% %                 % Retrieve points in a leaf node
% %                 leaf_cell_ptc_ids = vertcat(OQTR.cell_pts(bridge_component_cell_ids(j)).id);
% %                 leaf_cell_ptc_xyz = OQTR.pts(leaf_cell_ptc_ids,1:3);
% %                 
% %                 % Remove points to be assigned
% %                 mask = ismember(leaf_cell_ptc_ids,bridge_component_cell_ptc_ids);
% %                 cell_remain_ptc_ids = leaf_cell_ptc_ids(~mask);
% %                 cell_remain_ptc_xyz = leaf_cell_ptc_xyz(~mask,1:3);
% %                 
% %                 % Extract additional points of the slab above the bottom surface: 
% %                 mask = bridge_component_cell_surface_centroid(3) - 0.01 <= cell_remain_ptc_xyz(:,3);
% %                 bridge_slab_ptc_ids = cell_remain_ptc_ids(mask);
% %                 bridge_slab_ptc_xyz = cell_remain_ptc_xyz(mask,1:3);
% %                 bridge_sub_part_ptc_ids = cell_remain_ptc_ids(~mask);
% %                 if numel(bridge_slab_ptc_ids) > 10
% %                     bridge_slab.cell(count).id = bridge_component_cell_ids(j);
% %                     bridge_slab.cell(count).ptc_ids = bridge_slab_ptc_ids;
% %                     bridge_slab.cell(count).ptc_xyz = bridge_slab_ptc_xyz;
% %                     count = count + 1;
% %                     % Remove assigned points from the QOTR
% % %                     OQTR.cell_pts(bridge_component_cell_ids(j)).id = cell_remain_ptc_ids(~mask); % Remaining ptc_ids   
% %                 end
% %                 % Remove assigned points from the QOTR
% %                 OQTR.cell_pts(bridge_component_cell_ids(j)).id = bridge_sub_part_ptc_ids;
% %                 if isempty(bridge_sub_part_ptc_ids)
% %                     OQTR.cell_props(bridge_component_cell_ids(j)) = 0;
% %                 end
% %             end
% %         end
% %     end
% %     %
% %     % For the cells out of the bottom surface but would be footpath and
% %     % papapet
% %     % Retrieve the cells of the bottom surface of the slab
% %     bridge_bottom_slab_cell_ids = (arrayfun(@(x) vertcat(Region(x).cell.id),bridge_structure_region_ids,'UniformOutput',false));
% %     bridge_bottom_slab_cell_ids = cell2mat(bridge_bottom_slab_cell_ids);
% %     bridge_bottom_slab_cell_surface_centroid = (arrayfun(@(x) vertcat(Region(x).cell.surface_centroid),bridge_structure_region_ids,'UniformOutput',false));
% %     bridge_bottom_slab_cell_surface_centroid = cell2mat(bridge_bottom_slab_cell_surface_centroid);
% % 
% % %     non_bridge_component_cell_ids = setdiff(leaf_cell_ids,bridge_component_cell_ids);
% %     cur_horizon_level = length(Bridge_Component);
% %     % For the cells of the bottom surfaces
% %     for i=1:cur_horizon_level
% %         if (~isempty(cell2mat(strfind(Bridge_Component(i).Description, 'bridge footpath')))||...
% %             ~isempty(cell2mat(strfind(Bridge_Component(i).Description, 'Parapet'))))
% %             % Retrieve cell ids within the bridge component
% %             bridge_component_cell_ids = vertcat(Bridge_Component(i).cell.id);
% % 
% %             for j=1:numel(bridge_component_cell_ids)
% %                 % Check if the cells are not the bottom slab
% %                 % Retrieve points in a leaf node
% %                 leaf_cell_ptc_ids = vertcat(OQTR.cell_pts(bridge_component_cell_ids(j)).id);
% %                 leaf_cell_ptc_xyz = OQTR.pts(leaf_cell_ptc_ids,1:3);
% %                 if isempty(find(bridge_bottom_slab_cell_ids == bridge_component_cell_ids(j),1))
% %                     % Find the close cell in bridge desk
% %                     bridge_desk_cell_id = Closest_Cell(OQTR, bridge_component_cell_ids(j), bridge_bottom_slab_cell_ids);
% %                     mask = ismember(bridge_bottom_slab_cell_ids,bridge_desk_cell_id);
% %                     bridge_desk_cell_center_z = mean(bridge_bottom_slab_cell_surface_centroid(mask,3));
% %                 else
% %                     mask = ismember(bridge_bottom_slab_cell_ids,bridge_component_cell_ids(j));
% %                     bridge_desk_cell_center_z = mean(bridge_bottom_slab_cell_surface_centroid(mask,3));
% %                 end
% %                     
% %                 % Extract the points of the slab
% % %                     mask = bridge_bottom_desk_region_cent(3) <= leaf_cell_ptc_xyz(:,3);
% %                 mask = bridge_desk_cell_center_z - 0.05 <= leaf_cell_ptc_xyz(:,3);
% %                 bridge_super_structure_cell_ptc_ids = leaf_cell_ptc_ids(mask);
% %                 bridge_super_structure_cell_ptc_xyz = leaf_cell_ptc_xyz(mask,1:3);
% %                 bridge_sub_structure_cell_ptc_ids = leaf_cell_ptc_ids(~mask);
% % 
% %                 if numel(bridge_super_structure_cell_ptc_ids) > 10
% %                     bridge_slab.cell(count).id = bridge_component_cell_ids(j);
% %                     bridge_slab.cell(count).ptc_ids = bridge_super_structure_cell_ptc_ids;
% %                     bridge_slab.cell(count).ptc_xyz = bridge_super_structure_cell_ptc_xyz;
% %                     count = count + 1;
% %                     % Remove assigned points from the QOTR
% % %                     OQTR.cell_pts(bridge_component_cell_ids(j)).id = cell_remain_ptc_ids(~mask); % Remaining ptc_ids   
% %                 end
% %                 % Remove assigned points from the QOTR
% %                 OQTR.cell_pts(bridge_component_cell_ids(j)).id = bridge_sub_structure_cell_ptc_ids;
% %                 if numel(bridge_sub_structure_cell_ptc_ids) < 10
% %                     OQTR.cell_props(bridge_component_cell_ids(j)) = 0;
% %                 end
% % 
% %             end
% %         end
% %     end
% %     % Update bridge component: Slab Side
% %     side_slab_ptc_xyz = vertcat(bridge_slab.cell.ptc_xyz);
% %     min_boudingbox = min2DBoundingBox(side_slab_ptc_xyz(:,1:2)');
% %     cur_horizon_level = length(Bridge_Component);
% %     Bridge_Component(cur_horizon_level+1).Description = 'Side Slab';
% %     Bridge_Component(cur_horizon_level+1).cell = bridge_slab.cell;
% %     Bridge_Component(cur_horizon_level+1).features = [mean(side_slab_ptc_xyz, 1),min_boudingbox.area,min_boudingbox.long_edge,min_boudingbox.short_edge];
% %     clear bridge_slab side_slab_ptc_xyz min_boudingbox        
% 
% elseif horizontal_level == 3 % Ground
%     % Determine ground Points: ground cell if there is small number of point belows the cell
%     Ground = struct('cell',[]); 
%     count = 1;
%     for i = 1:length(Region)
%         
%         if Region_Features(i,12) >= 2 % ground region
%             
%             region_cell_ids = vertcat(Region(i).cell.id);
%             
%             for j = 1:numel(region_cell_ids)
%                 % retrieve points within the cell of the region
%                 region_cell_ptc_ids = Region(i).cell(j).ptc_ids;
%                 region_cell_surface_centroid = Region(i).cell(j).surface_centroid;
%                 % Retrieve points within the leaf node
%                 leaf_cell_ptc_ids = vertcat(OQTR.cell_pts(region_cell_ids(j)).id);
%                 non_footpath_cell_ptc_xyz = OQTR.pts(leaf_cell_ptc_ids,1:3);
% 
%                 % Remove points to be assigned
%                 mask = ismember(leaf_cell_ptc_ids,region_cell_ptc_ids);
%                 cell_remain_ptc_ids = leaf_cell_ptc_ids(~mask);
%                 cell_remain_ptc_xyz = non_footpath_cell_ptc_xyz(~mask,:);
%                 % Update group points
%                 
%                 % Insert the function to avod doplicate cell ids
%                 if ~isempty(Ground.cell)
%                     mask = vertcat(Ground.cell.id) == region_cell_ids(j);
%                     current_cell_order = find(mask == 1);
%                 else
%                     current_cell_order = [];
%                     
%                 end
%                 
%                 if isempty(current_cell_order)
%                     Ground.cell(count).id = Region(i).cell(j).id;
%                     Ground.cell(count).ptc_ids = Region(i).cell(j).ptc_ids;
%                     Ground.cell(count).ptc_xyz = Region(i).cell(j).ptc_xyz;
%                     Ground.cell(count).surface_centroid = Region(i).cell(j).surface_centroid;
%                     Ground.cell(count).surface_normls = Region(i).cell(j).surface_normls;
%                     Ground.cell(count).status = Region(i).cell(j).status;
%                     count = count + 1;
%                 else
%                     ground_cell_ptc_ids = Ground.cell(current_cell_order).ptc_ids;
%                     ground_cell_ptc_ids = union(ground_cell_ptc_ids, Region(i).cell(j).ptc_ids, 'rows');
%                     ground_cell_ptc_xyz = OQTR.pts(ground_cell_ptc_ids,1:3);
%                     [wptc_surface_centroid, wptc_surface_normals] = wPCA('ptc', ground_cell_ptc_xyz, 'ini_normal_vector', [0, 0, 1]);
%         
%                     % Update
% %                     Ground.cell(count).id = Region(i).cell(j).id;
%                     Ground.cell(current_cell_order).ptc_ids = ground_cell_ptc_ids;
%                     Ground.cell(current_cell_order).ptc_xyz = ground_cell_ptc_xyz;
%                     Ground.cell(current_cell_order).surface_centroid = wptc_surface_centroid;
%                     Ground.cell(current_cell_order).surface_normls = wptc_surface_normals;
%                     Ground.cell(current_cell_order).status = Region(i).cell(j).status;
% 
%                 end
% 
%                 % Update leaf node: Ground cell -> keep only points
%                 % above a local surface
%                 mask = cell_remain_ptc_xyz(:,3) >= region_cell_surface_centroid(3);
%                 OQTR.cell_pts(region_cell_ids(j)).id = cell_remain_ptc_ids(mask);
%                 if sum(mask) == 0
%                     OQTR.cell_props(region_cell_ids(j)) = 0;
%                 end     
%             end 
%         end
%     end
%  
% else %Ignore
% end



























%% Update OQTR
% For planar cell
planar_cell_ids = Planar_Cell.cell_ids;
for i = 1:numel(planar_cell_ids)
    mask = ismember(Planar_Cell.cell_ids, planar_cell_ids(i)); 
    planar_cell_ptc_ids =  Planar_Cell.ptc_ids(mask).id;
    OQTR.cell_pts(planar_cell_ids(i)).id = planar_cell_ptc_ids;
end

non_planar_cell_ids = nonPlanar_Cell.cell_ids;
for i = 1:numel(non_planar_cell_ids)
    mask = ismember(nonPlanar_Cell.cell_ids, non_planar_cell_ids(i)); 
    non_planar_cell_ptc_ids =  nonPlanar_Cell.ptc_ids(mask).id;
    OQTR.cell_pts(non_planar_cell_ids(i)).id = non_planar_cell_ptc_ids;
end

% Update properties of the cell
leaf_cell_ids = Node_Leaf(OQTR);
for i = 1:numel(leaf_cell_ids)
    if isempty(OQTR.cell_pts(leaf_cell_ids(i)).id) 
        i
        OQTR.cell_props(leaf_cell_ids(i)) = 0;
    end
end
%% Save data from segmentation
num_region = length(Region);
regions_pts = [];
for i = 1:num_region
   region_pts_xyz = vertcat(Region(i).cell.ptc_xyz(:,1:3));
   region_pts_xyz(:,4) = i;
   regions_pts = [regions_pts;region_pts_xyz];
   clear region_pts_xyz
end

outputFile = strcat(strtok(filename,'.'),'Segmentation.txt');
write_txt(outputFile, regions_pts(:,1:4))

%% Save data from segmentation and after refinement
num_region = length(Region);
regions_pts = [];
for i = 1:num_region
   region_pts_xyz = vertcat(Region(i).cell.ptc_xyz);
   region_add_ptc_xyz = unique(vertcat(Region(i).cell.add_ptc_xyz),'rows');
   region_pts_xyz = union(region_pts_xyz, region_add_ptc_xyz,'rows');
   region_pts_xyz(:,4:end) = [];
   region_pts_xyz(:, 4) = i;
   regions_pts = [regions_pts;region_pts_xyz];
   clear region_pts_xyz
end

outputFile = strcat(strtok(filename,'.'),'_Segmentation_Refinement.txt');
write_txt(outputFile, regions_pts(:,1:4))

%% Extract a point cloud of the bridge deck: large region and smallest deviation angle with nz
region_id = 1;
region_ptc_xyz = vertcat(Region(region_id).cell.ptc_xyz);
region_add_ptc_xyz = unique(vertcat(Region(region_id).cell.add_ptc_xyz),'rows');
bridge_deck_ptc_xyz = union(region_ptc_xyz, region_add_ptc_xyz,'rows');

% Write to text files
% outputFile = strcat(strtok(filename,'.'),'_CellSize_',num2str(THRESHOLD.cell_size*100),'_cm_backward.txt');
% write_txt(outputFile, region_ptc_xyz)
outputFile = strcat(strtok(filename,'.'),'Angle_',num2str(ceil(rad2deg(acos(THRESHOLD.max_angle)))),'_CellSize_',num2str(THRESHOLD.cell_size*1000),'_mm_',...
    'dev_height_', num2str(THRESHOLD.dev_height*1000), 'mm_bridgedeck.txt');
write_txt(outputFile, bridge_deck_ptc_xyz)
clear region_id region_ptc_xyz region_add_ptc_xyz outputFile
%%
% pts_id = vertcat(Leaf.ptc_id.id);
% pts_xyz = OQTR.Points(pts_id,:);
% write_txt('Leaf.txt', pts_xyz)



%% Write the result for evaluation
tic

outputFile = strcat(strtok(filename,'.'),'Angle_',num2str(ceil(rad2deg(acos(THRESHOLD.max_angle)))),'_CellSize_',num2str(THRESHOLD.cell_size*1000),'_mm_',...
    'dev_height_', num2str(THRESHOLD.dev_height*1000), 'mm_bridgedeck_Classification.txt');
mask = ismember(OQTR.Points, bridge_deck_ptc_xyz,  'rows');

write_txt(outputFile, [[OQTR.Points(mask,:),ones(sum(mask),1)];...
                       [OQTR.Points(~mask,:), zeros(sum(~mask),1)]])
toc
clear outputFile
%% Fitting the bridge deck
x = bridge_deck_ptc_xyz(:,1);
y = bridge_deck_ptc_xyz(:,2);
z = bridge_deck_ptc_xyz(:,3);
[fitobject, gof, output] = fit([x,y],z,'a + b*x + c*y + d*x.^2 + e*y.^2 + f*x.*y', 'StartPoint', [1, 1,1,1,1,1]);

z0 = fitobject(x, y);
results = z - z0;

outputFile = strcat(strtok(filename,'.'),'Angle_',num2str(ceil(rad2deg(acos(THRESHOLD.max_angle)))),'_CellSize_',num2str(THRESHOLD.cell_size*1000),'_mm_',...
    'dev_height_', num2str(THRESHOLD.dev_height*1000), 'mm_bridgedeck_diffElevation.txt');

write_txt(outputFile, [bridge_deck_ptc_xyz, results*1000]);
clear x y z z0 results
%%
% [fitobject, gof, output] = fit([x,y],z,'poly23');
% 
% z0 = fitobject(x, y);
% results = z - z0;
% 
% outputFile = strcat(strtok(filename,'.'),'_CellSize_',num2str(THRESHOLD.cell_size*100),'_cm_ele_deviation.txt');
% 
% write_txt(outputFile, [bridge_deck_ptc_xyz, results*1000]);
% clear x y z z0 results

%% Compute features of the points: Residuals, curvature and density
% http://www.academicjournals.org/article/article1380808355_Moazami%20et%20al.pdf
range_search_radius = 0.25;
[idx, dist]= rangesearch(bridge_deck_ptc_xyz, bridge_deck_ptc_xyz, range_search_radius);
ptc_features = zeros(size(bridge_deck_ptc_xyz,1),3);
%
for i = 1:size(bridge_deck_ptc_xyz,1)
    neighbour_ptc_id = idx{i};
    if numel(neighbour_ptc_id) >= THRESHOLD.min_num_pts
        neighbour_ptc_xyz = bridge_deck_ptc_xyz(neighbour_ptc_id,:);
        
        % Residual
        [~, eigvalue, ptc_residual]  = eigenspace(neighbour_ptc_xyz, 4);
        % Curvature
        ptc_curvature = eigvalue(1)/sum(eigvalue);
        % Density
        ptc_density = numel(dist{i})/(pi*max(dist{i})^2);
        ptc_features(i,:) = [ptc_density, ptc_residual, ptc_curvature];
    else
        ptc_features(i,:) = inf;
    end
end
clear range_search_radius idx dist i neighbour_ptc_id neighbour_ptc_xyz eigvalue ptc_residual ptc_curvature ptc_density 
%% Write data
mask = ~isinf(ptc_features(:,1));

outputFile = strcat(strtok(filename,'.'),'_CellSize_',num2str(THRESHOLD.cell_size*100),'_cm_bridge_deck_health.txt');

write_txt(outputFile, [bridge_deck_ptc_xyz(mask,:), ptc_features(mask,:)]);
clear mask
%% Write all regions
num_region = length(Region);
all_ptc = [];
for i = 1:num_region
    region_ptc_xyz = vertcat(Region(i).cell.ptc_xyz);
    region_add_ptc_xyz = unique(vertcat(Region(i).cell.add_ptc_xyz),'rows');
    region_all_ptc_xyz = union(region_ptc_xyz, region_add_ptc_xyz,'rows');
    region_all_ptc_xyz(:,4) = i;
    all_ptc = [all_ptc;region_all_ptc_xyz];
end
outputFile = strcat(strtok(filename,'.'),'_CellSize_',num2str(THRESHOLD.cell_size*100),'_cm.txt');

write_txt(outputFile, all_ptc);



%%







d = vertcat(Planar_Cell.ptc_xyz.id);
size(d)
size(unique(d,'rows'))

size(region_ptc_xyz)
size(unique(region_ptc_xyz,'rows'))

size(region_add_ptc_xyz)
size(unique(region_add_ptc_xyz,'rows'))



%% 
region_id = 2;
ptc = vertcat(Region(region_id).cell.ptc_xyz);
write_txt('deck_backward.txt', ptc)


%% Forward
ptc_2 = vertcat(Region(region_id).cell.add_ptc_xyz);
ptc_2 = unique(ptc_2, 'rows');
write_txt('deck_forward.txt', union(ptc, ptc_2,'rows'))



%%
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
hold all
ptc = vertcat(Planar_Cell.ptc_xyz.id);
Point3DPlot(ptc,7,'hsv', 1);



%%
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
hold all
% plot(zi, fz)
plot(fz, zi)
% plot(pks, zz)
%%

h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
hold all

plotAll(OQTR,'voxelId',region_cell_id,'color','r', 'fill',0,'typePlotPoint',1) 
%% Plot voxel
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
hold all
cell_id_1 = 99;
cell_id_2 = [98 100 101 104]';
plotAll(OQTR,'voxelId',cell_id_1,'color','r', 'fill',0) 
plotAll(OQTR,'voxelId',cell_id_2,'color','c','fill',0) 

mask = ismember(Planar_Cell.cell_id,cell_id_1,'rows');
ptc1 = Planar_Cell.ptc_xyz(mask).id;
plot3(ptc1(:,1),ptc1(:,2),ptc1(:,3),'r.', 'markersize',5);

mask = ismember(Planar_Cell.cell_id,cell_id_2,'rows');
ptc2 = vertcat(Planar_Cell.ptc_xyz(mask).id);
plot3(ptc2(:,1),ptc2(:,2),ptc2(:,3),'b.', 'markersize',5);
clear ptc1 ptc2 cell_id mask

%%
cell_id = 99;
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
hold all
plotAll(OQTR,'voxelId',cell_id_1,'color','r','typePlotPoint',1,'fill', 0) 

%% Plot region infor
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
hold all
region_no = 2;
mask = region_info(:,2) == region_no;
cell_ind = region_info(mask,1);
plotAll(OQTR,'voxelId',cell_ind,'color','r','typePlotPoint',1,'fill', 0) 

clear mask region_no cell_ind




%% Visualize the cell contains point clouds of multiple component
% nodeLeaf = leafNode(OQTR);

bandwidth = THRESHOLD.cell_size*THRESHOLD.max_slope;
filter_scale = 2;
%% 1010 Barie; 1020: Parapet + ground; 1005: parapet + pedesitrian; 1011: 1012; 1999
close all;
% i = 2560;
i=1999
cell_id = leaf_cell_ids(i);
cell_pts_ids = OQTR.voxelPoints(cell_id).id;
cell_pts_xyz = OQTR.Points(cell_pts_ids,:);

% Using smooth kernel to remove incorrect
[fz,zi,~] = ksdensity(cell_pts_xyz(:,3),'npoints',100,'bandwidth', bandwidth);

% Extract all peaks with a distance by the bandwidth
[pks,locs,~,~] = findpeaks(fz, 'MinPeakDistance', bandwidth);
zii = zi(locs);
[~, sort_id] = sort(zii,'descend');
zii = zii(sort_id);
pks = pks(sort_id);

for ii = 1: numel(pks)-1
    z_peak = zii(ii);
    if pks(ii) <= pks(ii+1)
        z_peak = zii(ii+1);
        break;
    end
    
end

% 
% %%
% 
% 
% % Filter the peaks larger than mean_peaks and extract the peak in the
% % highest elevation
% mean_peak = mean(fz);
% mask  = pks >= mean_peak;
% zii = zii(mask);
% [~, ind] = max(zii);
% z_peak = zii(ind);

% Extract the points within the peak
mask = (cell_pts_xyz(:,3) >= z_peak - filter_scale*bandwidth)&(cell_pts_xyz(:,3) <= z_peak + filter_scale*bandwidth);

node_leaf_pts_xyz_2 = cell_pts_xyz(mask,:);
    
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[0 0 800 800]);
% set(gca,'DataAspectRatio',[5 1 1]);

hold all
% plot(fz, zi, 'Linewidth', 2)
plot3(cell_pts_xyz(:,1), cell_pts_xyz(:,2), cell_pts_xyz(:,3),'.', 'markersize', 8)
x = (min(cell_pts_xyz(:,1)) + fz/15);
y = mean(cell_pts_xyz(:,2))*ones(numel(fz),1);
z = zi;
plot3(x, y, z, 'Linewidth', 2, 'color','b')

plot3(node_leaf_pts_xyz_2(:,1), node_leaf_pts_xyz_2(:,2), node_leaf_pts_xyz_2(:,3),'ro', 'markersize', 5)


grid on
view(30,15)
ax = gca;

ax.XAxis.TickLabelFormat = '%,.1f';
ax.YAxis.TickLabelFormat = '%,.1f';
ax.FontSize = 16;
ax.FontName = 'Time New Roman';
% write_txt('temp.txt', node_leaf_pts_xyz)
F = getframe(gcf);
imwrite(F.cdata,strcat('Multiple_point_in_Cell_',num2str(i),'_Cell_size_50cm.tif'),'Resolution',300); 

%% 1010 Barie; 1020: Parapet + ground; 1005: parapet + pedesitrian; 1011: 1012; 1999
leaf_cell_ids = leafNode(OQTR);
bandwidth = THRESHOLD.cell_size*THRESHOLD.max_slope;
filter_scale = 2;
%% 2640 2642
i = 2640
cell_id = leaf_cell_ids(i);
cell_pts_ids = OQTR.voxelPoints(cell_id).id;
cell_pts_xyz = OQTR.Points(cell_pts_ids,:);

% Using smooth kernel to remove incorrect
[fz,zi,~] = ksdensity(cell_pts_xyz(:,3),'npoints',100,'bandwidth', bandwidth);
write_txt('temp.txt', cell_pts_xyz)
%%
% Extract all peaks with a distance by the bandwidth
[pks,locs,~,~] = findpeaks(fz, 'MinPeakDistance', bandwidth);

zii = zi(locs);
[~, sort_id] = sort(zii,'descend');
zii = zii(sort_id);
pks = pks(sort_id);
if numel(pks) < 2
    z_peak = zii;
else
    z_peak = zii(1);
    p_peak = pks(1);
    for ii = 2: numel(pks)

        if (p_peak < 1.25*pks(ii))&&(z_peak - zii(ii) < 1)
            z_peak = zii(ii);
            p_peak = pks(ii);
        end
    end
end
%%
close all
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[0 0 800 800]);
% set(gca,'DataAspectRatio',[5 1 1]);

hold all
% plot(fz, zi, 'Linewidth', 2)
plot3(cell_pts_xyz(:,1), cell_pts_xyz(:,2), cell_pts_xyz(:,3),'.', 'markersize', 8)
x = (min(cell_pts_xyz(:,1)) + fz/15);
y = mean(cell_pts_xyz(:,2))*ones(numel(fz),1);
z = zi;
plot3(x, y, z, 'Linewidth', 2, 'color','b')
%%
mask = (cell_pts_xyz(:,3) >= z_peak - filter_scale*bandwidth)&(cell_pts_xyz(:,3) <= z_peak + filter_scale*bandwidth);

node_leaf_pts_xyz_2 = cell_pts_xyz(mask,:);

plot3(node_leaf_pts_xyz_2(:,1), node_leaf_pts_xyz_2(:,2), node_leaf_pts_xyz_2(:,3),'ro', 'markersize', 5)


grid on
view(30,15)
write_txt('temp.txt', node_leaf_pts_xyz_2)
%% Visualize backward algorithm
i = 1;
mask = region_info(:,2) == i;
region_cell_id = region_info(mask,1);

% Find cells on boundary of the region, which has less than 8 neibour
% region
neighbour_cell = struct('id',[]);
flag = zeros(numel(region_cell_id),2); %[cell_id, 1: bound or 0: interior]
flag(:,1) = region_cell_id;

for j = 1: numel(region_cell_id)
    neighbour_cell_id = queryVoxelNeighbour(TREE, region_cell_id(j), region_cell_id);
    neighbour_cell(j).id = neighbour_cell_id;
    if numel(neighbour_cell_id) < 8 %Bound cell

        flag(j,2) = 1; %Bound cell
    end
end

%% 50 51: forward
j = 118;
% Current cell in a region
region_current_cell_id = region_cell_id(j);
mask = ismember(LeafNode.cell_id, region_current_cell_id);
current_cell_ptc_xyz = vertcat(LeafNode.ptc_xyz(mask).id);
% size(current_cell_ptc_xyz)
% write_txt('temp.txt', current_cell_ptc_xyz)

bound_cell_id = flag(logical(flag(:,2)),1);
neighbour_cell_id = neighbour_cell(j).id;
mask = ismember(neighbour_cell_id, bound_cell_id,'rows');

interior_cell_id = neighbour_cell_id(~mask);

exterior_cell_id = neighbour_cell_id(mask);

all_interior_cell_id = interior_cell_id;
for k = 1:numel(interior_cell_id)
    neighbour_interior_cell_id = queryVoxelNeighbour(TREE, interior_cell_id(k), region_cell_id);
    all_interior_cell_id = union(all_interior_cell_id, neighbour_interior_cell_id, 'rows');
end
% Remove bound_cell again
mask = ismember(all_interior_cell_id, bound_cell_id,'rows');

all_interior_cell_id = all_interior_cell_id(~mask);

% Searching neighbour cells not in the region
out_neighbour_cell_id = queryVoxelNeighbour(TREE, region_current_cell_id, leaf_cell_ids);
        
% Remove neighbour_bound_cell_id belongs to the current region
mask = ismember(out_neighbour_cell_id, region_cell_id, 'rows');
out_neighbour_cell_id = out_neighbour_cell_id(~mask);

out_neighbour_bound_cell_ptc_id = vertcat(TREE.voxelPoints(out_neighbour_cell_id).id);
out_neighbour_bound_cell_ptc_xyz = TREE.Points(out_neighbour_bound_cell_ptc_id,:);



%% Plot current cell, interior and exterior cells
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[0 0 800 800]);
hold all
% Plot points of interest 
voxelNeighbourId = windowQueryVoxelNeighbour(OQTR, region_current_cell_id, region_cell_id, 1.5);
mask = ismember(voxelNeighbourId,[2696; 2694;2636;2637; 2634; 2635;2631;2632;2633],'rows');
voxelNeighbourId = voxelNeighbourId(~mask);
mask = ismember(LeafNode.cell_id,voxelNeighbourId);

interest_ptc_id = vertcat(LeafNode.ptc_id(mask).id);
interest_ptc_xyz = vertcat(OQTR.Points(interest_ptc_id,1:3));
plot3(interest_ptc_xyz(:,1), interest_ptc_xyz(:,2), interest_ptc_xyz(:,3),'k.', 'markersize', 8)

% plot point out of the region
out_neighbour_cell_id_2 = windowQueryVoxelNeighbour(OQTR, out_neighbour_cell_id(1), leaf_cell_ids, 1);
mask = ismember(out_neighbour_cell_id_2, region_cell_id, 'rows');
out_neighbour_cell_id_2 = out_neighbour_cell_id_2(~mask);

mask = ismember(LeafNode.cell_id,out_neighbour_cell_id_2);

out_ptc_id = vertcat(LeafNode.ptc_id(mask).id);
out_ptc_xyz = vertcat(OQTR.Points(out_ptc_id,1:3));

plot3(out_ptc_xyz(:,1), out_ptc_xyz(:,2), out_ptc_xyz(:,3),'b.', 'markersize', 8)

% Plot the current cell
temp_tree = OQTR;
temp_tree.voxelBoundaries(:,3) = temp_tree.voxelBoundaries(:,6)-0.1;
plotAll(temp_tree,'voxelId',region_current_cell_id,'color','r','typePlotPoint',0,'fill', 0)

% Interior cell 
plotAll(temp_tree,'voxelId',interior_cell_id,'color','m','typePlotPoint',0,'fill', 0)

% Plot exterior cells
plotAll(temp_tree,'voxelId',exterior_cell_id,'color','c','typePlotPoint',0,'fill', 0)

% Plot out region cell
for kk = 1:numel(out_neighbour_cell_id)
    mask = ismember(LeafNode.cell_id,out_neighbour_cell_id(kk));
    out_ptc_id = vertcat(LeafNode.ptc_id(mask).id);
    out_ptc_xyz = vertcat(OQTR.Points(out_ptc_id,1:3));
    temp_tree.voxelBoundaries(out_neighbour_cell_id(kk),3) = min(out_ptc_xyz(:,3));
    temp_tree.voxelBoundaries(out_neighbour_cell_id(kk),6) = max(out_ptc_xyz(:,3));
end


plotAll(temp_tree,'voxelId',out_neighbour_cell_id,'color','g','typePlotPoint',0,'fill', 0)

% plotAll(OQTR,'voxelId',region_current_cell_id,'color','k','typePlotPoint',0,'fill', 0)
% plotAll(OQTR,'voxelId',interior_cell_id,'color','r','typePlotPoint',0,'fill', 0)
% 

%%
close all;
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[0 0 800 800]);
set(gca,'DataAspectRatio',[1 .75 1]);
hold all
% Plot points of interest 
voxelNeighbourId = windowQueryVoxelNeighbour(OQTR, region_current_cell_id, region_cell_id, 1.5);
mask = ismember(voxelNeighbourId,[2696; 2694;2636;2637; 2634; 2635;2631;2632;2633],'rows');
voxelNeighbourId = voxelNeighbourId(~mask);
mask = ismember(LeafNode.cell_id,voxelNeighbourId);

interest_ptc_id = vertcat(LeafNode.ptc_id(mask).id);
interest_ptc_xyz = vertcat(OQTR.Points(interest_ptc_id,1:3));
h1 = plot3(interest_ptc_xyz(:,1), interest_ptc_xyz(:,2), interest_ptc_xyz(:,3),'.', 'markersize', 8,'color',[0.5,0.5,0.5]);

% plot point out of the region
out_neighbour_cell_id_2 = windowQueryVoxelNeighbour(OQTR, out_neighbour_cell_id(1), leaf_cell_ids, 1);
mask = ismember(out_neighbour_cell_id_2, region_cell_id, 'rows');
out_neighbour_cell_id_2 = out_neighbour_cell_id_2(~mask);

mask = ismember(LeafNode.cell_id,out_neighbour_cell_id_2);

out_ptc_id = vertcat(LeafNode.ptc_id(mask).id);
out_ptc_xyz = vertcat(OQTR.Points(out_ptc_id,1:3));

h2 = plot3(out_ptc_xyz(:,1), out_ptc_xyz(:,2), out_ptc_xyz(:,3),'.', 'markersize', 8, 'color','b');

view(45,45)
ax = gca;

ax.XAxis.TickLabelFormat = '%,.1f';
ax.YAxis.TickLabelFormat = '%,.1f';
ax.FontSize = 18;
ax.FontName = 'Time New Roman';
% write_txt('temp.txt', node_leaf_pts_xyz)
lgd = legend([h1, h2],'A point cloud of R_i','A point cloud out of R_i');
%
lgd.Location = 'best';
lgd.Orientation = 'horizontal';
lgd.Box = 'off'; 
lgd.FontSize = 18;
lgd.Position = [0.4,0.5, 0.2, 0.6];

F = getframe(gcf);
imwrite(F.cdata,strcat('Point_interest.tif'),'Resolution',300); 

% ch = findobj(get(lgd,'children'), 'type', 'line'); %// children of legend of type line
% set(ch, 'Markersize', 24)
%%
close all;
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[0 0 800 800]);
set(gca,'DataAspectRatio',[1 .75 1]);
hold all

% points in a current 
mask = ismember(LeafNode.cell_id,region_current_cell_id);
current_cell_ptc_id = vertcat(LeafNode.ptc_id(mask).id);
current_cell_ptc_xyz = vertcat(OQTR.Points(current_cell_ptc_id,1:3));


% Plot points of interest 
voxelNeighbourId = windowQueryVoxelNeighbour(OQTR, region_current_cell_id, region_cell_id, 1.5);
mask = ismember(voxelNeighbourId,[2696; 2694;2636;2637; 2634; 2635;2631;2632;2633],'rows');
voxelNeighbourId = voxelNeighbourId(~mask);
mask = ismember(LeafNode.cell_id,voxelNeighbourId);

interest_ptc_id = vertcat(LeafNode.ptc_id(mask).id);
interest_ptc_xyz = vertcat(OQTR.Points(interest_ptc_id,1:3));
interest_ptc_xyz = setdiff(interest_ptc_xyz ,current_cell_ptc_xyz,'rows');

h1 = plot3(interest_ptc_xyz(:,1), interest_ptc_xyz(:,2), interest_ptc_xyz(:,3),'.', 'markersize', 8,'color',[0.5,0.5,0.5]);

% plot point out of the region
out_neighbour_cell_id_2 = windowQueryVoxelNeighbour(OQTR, out_neighbour_cell_id(1), leaf_cell_ids, 1);
mask = ismember(out_neighbour_cell_id_2, region_cell_id, 'rows');
out_neighbour_cell_id_2 = out_neighbour_cell_id_2(~mask);

mask = ismember(LeafNode.cell_id,out_neighbour_cell_id_2);

out_ptc_id = vertcat(LeafNode.ptc_id(mask).id);
out_ptc_xyz = vertcat(OQTR.Points(out_ptc_id,1:3));

h2 = plot3(out_ptc_xyz(:,1), out_ptc_xyz(:,2), out_ptc_xyz(:,3),'b.', 'markersize', 8);


h3 = plot3(current_cell_ptc_xyz(:,1), current_cell_ptc_xyz(:,2), current_cell_ptc_xyz(:,3),'b.', 'markersize', 7,...
    'color',[255/255,102/255,255/255]);


% Plot the current cell
temp_tree = OQTR;
temp_tree.voxelBoundaries(:,3) = temp_tree.voxelBoundaries(:,6)-0.1;
plotAll(temp_tree,'voxelId',region_current_cell_id,'color','r','typePlotPoint',0,'fill', 0);

% Interior cell 
plotAll(temp_tree,'voxelId',interior_cell_id,'color','m','typePlotPoint',0,'fill', 0);

% Plot exterior cells
plotAll(temp_tree,'voxelId',exterior_cell_id,'color','c','typePlotPoint',0,'fill', 0);

% Plot out region cell
for kk = 1:numel(out_neighbour_cell_id)
    mask = ismember(LeafNode.cell_id,out_neighbour_cell_id(kk));
    out_ptc_id = vertcat(LeafNode.ptc_id(mask).id);
    out_ptc_xyz = vertcat(OQTR.Points(out_ptc_id,1:3));
    temp_tree.voxelBoundaries(out_neighbour_cell_id(kk),3) = min(out_ptc_xyz(:,3));
    temp_tree.voxelBoundaries(out_neighbour_cell_id(kk),6) = max(out_ptc_xyz(:,3));
end


plotAll(temp_tree,'voxelId',out_neighbour_cell_id,'color','g','typePlotPoint',0,'fill', 0)


view(45,45)
ax = gca;

ax.XAxis.TickLabelFormat = '%,.1f';
ax.YAxis.TickLabelFormat = '%,.1f';
ax.FontSize = 18;
ax.FontName = 'Time New Roman';
% write_txt('temp.txt', node_leaf_pts_xyz)
%%
% legend_str{1} = '';
% legend_str{2} = '';
% legend_str{3} = 'Current cell C_i';
% legend_str{4} = 'Interior cell C_k';
% legend_str{5} = 'Boundary cell C_l';
% legend_str{6} = 'Cell C_m out of R_i';
lgd = legend('','',' A point cloud of C_i',' Current cell C_i',' Interior cell C_k', ' Boundary cell C_l',...
    ' Cell C_m out of R_i');
% lgd = columnlegend(3, legend_str)
%
lgd.Location = 'northeast';

lgd.Orientation = 'vertical';
% lgd.NumColumnsMode
% lgd.NumColumns = 2;
lgd.Box = 'off'; 
lgd.FontSize = 18;
lgd.Position = [0.7,0.67, 0.2, 0.2];
%
F = getframe(gcf);
imwrite(F.cdata,strcat('Point_interest_cell.tif'),'Resolution',300); 

% ch = findobj(get(lgd,'children'), 'type', 'line'); %// children of legend of type line
% set(ch, 'Markersize', 24)%
%% Plot the point to be remove
close all;
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[0 0 650 800]);
set(gca,'DataAspectRatio',[1 .75 1]);
hold all

% Remain points after backward 
mask = ismember(vertcat(Region(1).cell.id),region_current_cell_id);

remain_ptc_xyz = Region(1).cell(mask).ptc_xyz(:,1:3);

h2 = plot3(remain_ptc_xyz(:,1), remain_ptc_xyz(:,2), remain_ptc_xyz(:,3),'.', 'markersize', 10,...
    'color','k');

% points in a current 
mask = ismember(LeafNode.cell_id,region_current_cell_id);
current_cell_ptc_id = vertcat(LeafNode.ptc_id(mask).id);
current_cell_ptc_xyz = vertcat(OQTR.Points(current_cell_ptc_id,1:3));

% x_lim = [min(current_cell_ptc_xyz(:,1)), max(current_cell_ptc_xyz(:,1))];
% y_lim = [min(current_cell_ptc_xyz(:,2)), max(current_cell_ptc_xyz(:,2))];
% z_lim = [min(current_cell_ptc_xyz(:,3)), max(current_cell_ptc_xyz(:,3))];

current_cell_ptc_xyz = setdiff(current_cell_ptc_xyz,remain_ptc_xyz,'rows'); 
h1 = plot3(current_cell_ptc_xyz(:,1), current_cell_ptc_xyz(:,2), current_cell_ptc_xyz(:,3),'.', 'markersize', 15,...
    'color',[255/255,102/255,255/255]);

view(45,45)
ax = gca;

ax.XAxis.TickLabelFormat = '%,.1f';
ax.YAxis.TickLabelFormat = '%,.1f';
ax.FontSize = 14;
ax.FontName = 'Time New Roman';

% ax.XLim = x_lim;
y_lim = ax.YLim;
ax.YLim = [y_lim(1)+0.1, y_lim(2)];
% xlabel('x coordinate')
% ylabel('y coordinate')
% zlabel('z coordinate')
% ax.XLim = z_lim;
%
lgd = legend(' A point cloud removed out of C_i',' A remaining point cloud of C_i');
% lgd = columnlegend(3, legend_str)
%
lgd.Location = 'northeast';

lgd.Orientation = 'horizontal';
% lgd.NumColumnsMode
% lgd.NumColumns = 2;
lgd.Box = 'off'; 
lgd.FontSize = 14;
lgd.Position = [0.4,0.625, 0.2, 0.2];
%
F = getframe(gcf);
imwrite(F.cdata,strcat('Point_backward.tif'),'Resolution',300); 

%% Forward

close all;
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[0 0 600 800]);
set(gca,'DataAspectRatio',[1 .75 1]);
hold all



% points in a current 
% mask = ismember(LeafNode.cell_id,region_current_cell_id);
% current_cell_ptc_id = vertcat(LeafNode.ptc_id(mask).id);
% current_cell_ptc_xyz = vertcat(OQTR.Points(current_cell_ptc_id,1:3));
% 
% current_cell_ptc_xyz = setdiff(current_cell_ptc_xyz,remain_ptc_xyz,'rows'); 
mask = ismember(vertcat(Region(1).cell.id),region_current_cell_id);

remain_ptc_xyz = Region(1).cell(mask).ptc_xyz(:,1:3);

h2 = plot3(remain_ptc_xyz(:,1), remain_ptc_xyz(:,2), remain_ptc_xyz(:,3),'.', 'markersize', 10,...
    'color',[255/255,102/255,255/255]);

% Remain points after backward 
mask = ismember(vertcat(Region(1).cell.id),region_current_cell_id);

add_ptc_xyz = Region(1).cell(mask).add_ptc_xyz(:,1:3);

h1 = plot3(add_ptc_xyz(:,1), add_ptc_xyz(:,2), add_ptc_xyz(:,3),'.', 'markersize', 13,...
    'color','k');

view(45,45)
ax = gca;

ax.XAxis.TickLabelFormat = '%,.1f';
ax.YAxis.TickLabelFormat = '%,.1f';
ax.FontSize = 14;
ax.FontName = 'Time New Roman';

% ax.XLim = x_lim;
y_lim = ax.YLim;
ax.YLim = [y_lim(1)+0.1, y_lim(2)];
% xlabel('x coordinate')
% ylabel('y coordinate')
% zlabel('z coordinate')
% ax.XLim = z_lim;
%%
lgd = legend(' A remaining point cloud of C_i',' A point cloud merged to R_i');
% lgd = columnlegend(3, legend_str)
%
lgd.Location = 'northwest';

lgd.Orientation = 'horizontal';
% lgd.NumColumnsMode
% lgd.NumColumns = 2;
lgd.Box = 'off'; 
lgd.FontSize = 14;
lgd.Position = [0.375,0.65, 0.2, 0.075];
%
F = getframe(gcf);
imwrite(F.cdata,strcat('Point_forward.tif'),'Resolution',300); 

