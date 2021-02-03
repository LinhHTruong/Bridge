%% Read file
clc; close all; clear all;
% curdir = 'd:\Point Cloud Processing\Point Cloud Data\Germany Concrete Bridge\';
curdir = 'c:\Users\ltruonghong\TUDelft\Point Cloud Processing\Data\UCD\Germany Concrete Bridge\';
% filename = 'Concrete_Bridge_10mm.txt';
filename = 'Sub_03_Concrete_Bridge_10mm.txt';

delimiter = {' ' '\b' '\t' ',' ';' }; % all possible delimiter 
fid = fopen(strcat(curdir,filename));
line = fgetl(fid);
firstLine = textscan(line, '%f', 'Delimiter', delimiter, 'MultipleDelimsAsOne', 1);
nocolumn = length(firstLine{1});
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
% Structural thresholds estblishing based on structural components
struct_threshold.bridge_min_length = 2;%10;
struct_threshold.footpath_min_width = 0.8;
struct_threshold.footpath_length_ratio = 0.25;
struct_threshold.parapet_min_height = 0.5;
struct_threshold.curb_max_height = 0.3; % Maximum height of the road curb, which mean the different height between the road surface and footpath is less than this threshold
struct_threshold.column_min_height = 0.5;
struct_threshold.adjoin_components_min_edge_length = 2.0; % The minimum length of edges sharing by two adjoin or parallel surfaces
struct_threshold.surface_error = 0.005;
struct_threshold.parapet_height_tol = 0.2;
struct_threshold.min_no_cells = 20;
struct_threshold.bandwidth = 0.1;
struct_threshold.search_scale = 2; % Window size to search neighbour cells
struct_threshold.section_length = 1.0;  % Minimum length of a cross-section
struct_threshold.section_height = 0.4; % Minimum height of a cross-section
struct_threshold.section_width = 0.4;  % Minimum width of a cross-section
struct_threshold.traffic_lane_width = 3.0; % Minimum width of a traffic lane
struct_threshold.section_interval = 0.2; % An interval 
struct_threshold.inclined_member = 45; % An inclinded angle of the structure in vertical direction 


% Threshold for filtering data 
sampling_step = 0.005;      % m, Real sampling step of the data
cell_size = 0.3;            % m, Input cell size and will be adjust with real cell size
min_num_cell = 3;           % The minimum number of the cells in the segment
min_num_pts = 10;           % The minimum number of the cells to estimate the plane
max_angle = 5.0;            % degree, Maximum deviation angle between two surfaces (local or global)
max_residual = 0.005;       % m, Maximum residual of the (local) surface for segmentation
max_distance = 0.01;        % m, Maximum distance of the (local) surface for segmentation
min_dist_2_planes = 0.2;    % m, A maximum distance between two parallel surfaces/planes
hor_plane_max_angle = 60;   % degree, The maximum angle between the horizontal plane and vertical
voxel_size = 0.05;          % m, Input voxel size

THRESHOLD = Sematic_thresholds('sampling_step',sampling_step,'cell_size',cell_size, 'min_num_cell',min_num_cell,'min_num_pts', min_num_pts,...
                                'max_angle', max_angle,'max_residual', max_residual,'max_distance', max_distance,...
                                'min_dist_2_planes', min_dist_2_planes,'hor_plane_max_angle',hor_plane_max_angle,'voxel_size',voxel_size);
clear sampling_step cell_size min_num_cell min_num_pts max_angle max_residual max_distance min_dist_2_planes hor_plane_max_angle voxel_size   

% Generate cell grids
tic
OQTR = OctQuadtree(PTC.xyz,'max_size', THRESHOLD.cell_size);
fprintf('Running time of generating quadtree: %.2f seconds \n', toc);
% clear PTC

% Extract the cells on a leaf node containing the horizontal plane
leaf_cell_ids = Node_Leaf(OQTR);
leaf_cell_ids = leaf_cell_ids(OQTR.cell_props(leaf_cell_ids) == 1);
% smallest_cell_size = (OQTR.cell_bounds(leaf_cell_ids,[4,5]) - OQTR.cell_bounds(leaf_cell_ids,[1,2]));
% THRESHOLD.cell_size = max(smallest_cell_size(:));

%%
[Cells_Plane, Cells_nonPlane] = cell_plane_extraction(OQTR, leaf_cell_ids, THRESHOLD);
clear smallest_cell_size

% Bridge structure#
    
time_record = true;
% Create data structure stored point clouds of the bridge components
SuperStructure = struct('Component',[], 'Linked_Component',{}, 'Description', {});

superstructure_name = {'top_surface', 'bottom_surface', 'intermediate_surface'};
superstructure_comp_name = {};
for bridge_tree_level = 1:2
    
    % Extract the top surface of the bridges
    [Cells_Plane, Cells_nonPlane, Cells_Region] = bridge_superstructure_surface_extraction(OQTR, Cells_Plane, Cells_nonPlane, THRESHOLD, 'top', true);

    % Update bridge tree structure

    if strcmp(superstructure_name{bridge_tree_level},'top_surface')
        % Road, footpath
        [SuperStructure, Cells_Region, un_connect_segs_ids, Road_Surface] = bridge_top_superstructure(OQTR, Cells_Region, SuperStructure, struct_threshold, THRESHOLD, time_record);
        
        %Road curb
        SuperStructure = extract_road_curb(OQTR, SuperStructure, THRESHOLD);
        
        %Parapet
        SuperStructure = extract_parapets(OQTR, SuperStructure, Road_Surface, struct_threshold, THRESHOLD);
        
    elseif strcmp(superstructure_name{bridge_tree_level},'bottom_surface')
        
        % Bottom surfaces and its connections
        [SuperStructure, Cells_Region, un_connect_segs_ids] = bridge_bottom_superstructure(OQTR, Cells_Region, SuperStructure, struct_threshold, THRESHOLD, time_record);

    end
    % Update cells
    Cells_Plane = update_cells(OQTR, Cells_Plane, Cells_Region, un_connect_segs_ids, SuperStructure, superstructure_comp_name, THRESHOLD, time_record);
%     Cells_Plane = update_cells(OQTR, Cells_Plane, Cells_Region, SuperStructure, superstructure_comp_name, un_connect_segs_ids, THRESHOLD, time_record);
%     
    % Update the name of the components
    superstructure_comp_name = [superstructure_comp_name;cat(1, SuperStructure.Description{:})];
end

%% Intermediate surface
[SuperStructure, Active_Pts] = bridge_intermediate_superstructure(OQTR, SuperStructure, struct_threshold, THRESHOLD, time_record);


%% Extract the clusters of abutments and pier
subOQTR = substructure_extraction (OQTR, Active_Pts, Road_Surface, struct_threshold, THRESHOLD);


%% Abutment extraction
Bridge_Abutment = bridge_abutment_extraction(OQTR, subOQTR, Road_Surface, struct_threshold, THRESHOLD, time_record);

%% Pier extraction
tic
Bridge_Pier = bridge_pier_extraction(OQTR, subOQTR, struct_threshold, THRESHOLD, time_record);
toc

%% Write the results
% Supstructure
num_comps = length(SuperStructure.Component);
for i = 1:num_comps
    % Retrieve the ptc_ids of each component
    comp_ptc_ids = vertcat(SuperStructure.Component(i).cell.ptc_ids);
    comp_ptc_xyz = OQTR.pts(comp_ptc_ids,1:3);
    comp_name = SuperStructure.Description{i};
    
    % Write txt
    out_file_name = strcat(comp_name, '.txt');
    write_txt(out_file_name, comp_ptc_xyz)
end
clear num_comps comp_ptc_ids comp_ptc_xyz
%% Write Abutment
num_abutments = length(Bridge_Abutment);
for i = 1:num_abutments
    % Retrieve the points in each surface
    abutment_name = Bridge_Abutment(i).Description;
    abutment_ptc_ids = vertcat(Bridge_Abutment(i).plane.ptc_ids);
    abutment_ptc_xyz = OQTR.pts(abutment_ptc_ids,1:3);
    
    % get plane id
    num_planes = length(Bridge_Abutment(i).plane);
    ptc_plane_id = inf(numel(abutment_ptc_ids),1);
    for j = 1:numel(num_planes)
        % Retrieve the ptc_ids of the plane
        plane_ptc_ids = Bridge_Abutment(i).plane(j).ptc_ids;
        mask = ismember(abutment_ptc_ids, plane_ptc_ids);
        ptc_plane_id(mask) = j;
    end
    % Write txt
    out_file_name = strcat(abutment_name, '.txt');
    write_txt(out_file_name, [abutment_ptc_xyz, ptc_plane_id])
end

%% Write Pier
num_piers = length(Bridge_Pier);
for i = 1:num_piers
    % Retrieve the points in each surface
    pier_name = Bridge_Pier(i).Description;
    pier_ptc_ids = vertcat(Bridge_Pier(i).plane.ptc_ids);
    pier_ptc_xyz = OQTR.pts(pier_ptc_ids,1:3);
    
    % get plane id
    num_planes = length(Bridge_Pier(i).plane);
    ptc_plane_id = inf(numel(pier_ptc_ids),1);
    for j = 1:numel(num_planes)
        % Retrieve the ptc_ids of the plane
        plane_ptc_ids = Bridge_Pier(i).plane(j).ptc_ids;
        mask = ismember(pier_ptc_ids, plane_ptc_ids);
        ptc_plane_id(mask) = j;
    end
    % Write txt
    out_file_name = strcat(pier_name, '.txt');
    write_txt(out_file_name, [pier_ptc_xyz, ptc_plane_id])
end
