
% Write point cloud of the 
% our_dir = 'D:\Point Cloud Processing\Matlab Codes\Bridge decomposition\Version 4.2\Version 42_';
out_dir = 'C:\Users\ltruonghong\TUDelft\Point Cloud Processing\Matlab codes\Bridge segmentation\Version 4.2\';
%% Write points with cell ids
tic
flag = zeros(size(OQTR.pts,1), 1);
for i=1:numel(leaf_cell_ids)
    leaf_cell_ptc_id = OQTR.cell_pts(leaf_cell_ids(i)).id;
    flag(leaf_cell_ptc_id,1) = leaf_cell_ids(i);
end
outputFile = strcat(out_dir, strtok(filename,'.'),'_LeafNode_Cell_ids.txt');
write_txt(outputFile, [OQTR.pts, flag(:,1)])
fprintf('Running time of writing remaining point cloud for next iteration: %.2f seconds \n', toc);
clear leaf_cells_ptc
%% Write results of plane
tic
flag = inf(size(OQTR.pts,1), 3);
flag(:,1) = 1: size(OQTR.pts,1);
for i = 1:length(Cells_Plane.cell_ids)
    
    if Cells_Plane.cell_ids(i,3) == 1

        cell_ptc_ids = Cells_Plane.peak_info(i).ptc_ids;
        cell_peak_id = Cells_Plane.cell_ids(i,:);
        flag(cell_ptc_ids,2) = cell_peak_id(1);
        flag(cell_ptc_ids,3) = cell_peak_id(2);
    end
end
mask = any(isinf(flag),2);
flag(mask,:) = [];

outputFile = strcat(our_dir, strtok(filename,'.'),'_Cells_Plane_Peaks_3.txt');
write_txt(outputFile, [OQTR.pts(flag(:,1),1:3), flag(:,2:3)])
fprintf('Running time of writing results of point extraction: %.2f seconds \n', toc);
clear cell_peak_id flag mask

%% Write results of plane
tic
flag = inf(size(OQTR.pts,1), 3);
flag(:,1) = 1: size(OQTR.pts,1);
for i = 1:length(cells_plane.cell_ids)
    
    if cells_plane.cell_ids(i,3) == 1

        cell_ptc_ids = cells_plane.peak_info(i).ptc_ids;
        cell_peak_id = cells_plane.cell_ids(i,:);
        flag(cell_ptc_ids,2) = cell_peak_id(1);
        flag(cell_ptc_ids,3) = cell_peak_id(2);
    end
end
mask = any(isinf(flag),2);
flag(mask,:) = [];

outputFile = strcat(our_dir, strtok(filename,'.'),'_Cells_Plane_Peaks_1.txt');
write_txt(outputFile, [OQTR.pts(flag(:,1),1:3), flag(:,2:3)])
fprintf('Running time of writing results of point extraction: %.2f seconds \n', toc);
clear cell_peak_id flag mask
%% Write data points of non_planar 
all_component_ptc_xyz = [];
for i = 1:length(nonPlane_Cells.cell_ids)
    component_ptc_ids = nonPlane_Cells.ptc_ids(i).id;
    component_ptc_xyz = OQTR.pts(component_ptc_ids,1:3);
    component_ptc_xyz(:,4) = nonPlane_Cells.cell_ids(i);
    all_component_ptc_xyz = [all_component_ptc_xyz;component_ptc_xyz];
end
outputFile = strcat(our_dir, strtok(filename,'.'),'_nonPlanar_Cell_1.txt');
write_txt(outputFile, all_component_ptc_xyz)
fprintf('Running time of writing results of point extraction: %.2f seconds \n', toc);
clear all_component_ptc_xyz component_ptc_ids component_ptc_xyz i
%% Center of all cell peak for segmentation
outputFile = strcat(our_dir, strtok(filename,'.'),'_cell_peak_center_1.txt');
write_txt(outputFile, cells_peaks_features(:,[3,4,5,end]))
clear no_region region_ptc_xyz mask current_cell_ptc_xyz outputFile

%% Write results of segmentation 
tic
horizontal_level = 1;
no_region = max(Region_Info(:,end));
all_region_ptc_xyz = [];
for i = 1:no_region
    mask = Region_Info(:,end) == i;
    region_cell_ids = Region_Info(mask,[1,2]);
    
    region_cell_ptc_xyz = [];
    for j = 1:size(region_cell_ids,1)
        cell_id = region_cell_ids(j,:);
        mask = ismember(Plane_Cells.cell_ids(:,[1,2]), cell_id, 'rows');
        cell_peak_ptc_ids = Plane_Cells.peak_info(mask).ptc_ids;
        
        % Extract points
        cell_peak_ptc_xyz = OQTR.pts(cell_peak_ptc_ids,1:3);
        
        cell_peak_ptc_xyz(:,4) = cell_id(1);
        region_cell_ptc_xyz = [region_cell_ptc_xyz;cell_peak_ptc_xyz];
    end
    region_cell_ptc_xyz(:,5) = i;
    all_region_ptc_xyz = [all_region_ptc_xyz;region_cell_ptc_xyz];

end

outputFile = strcat(our_dir, strtok(filename,'.'),'_Segmentation_1.txt');
write_txt(outputFile, all_region_ptc_xyz)
clear region_cell_ptc_xyz no_region all_region_ptc_xyz i region_cell_ids cell_id cell_peak_ptc_ids cell_peak_ptc_xyz outputFile

fprintf('Running time of writing results of point extraction: %.2f seconds \n', toc);
% %% Segmenttation with segment ids
% tic
% no_region = max(Region_Info(:,end));
% all_region_ptc_xyz = [];
% for i = 1:no_region
%     mask = Region_Info(:,end) == i;
%     region_cell_ids = Region_Info(mask,[1,2]);
%     mask = ismember(Plane_Cells.cell_ids(:,[1,2]), region_cell_ids, 'rows');
%     region_cell_peak_ptc_ids = vertcat(Plane_Cells.peak_info(mask).ptc_ids);
%     region_cell_peak_ptc_xyz = OQTR.pts(region_cell_peak_ptc_ids,1:3);
%     region_cell_peak_ptc_xyz(:,4) = i;
%     
%     all_region_ptc_xyz = [all_region_ptc_xyz;region_cell_peak_ptc_xyz];
% 
% end
% 
% outputFile = strcat(our_dir, strtok(filename,'.'),'_Segmentation_2.txt');
% write_txt(outputFile, all_region_ptc_xyz)
% clear region_cell_ptc_xyz no_region all_region_ptc_xyz i region_cell_ids cell_id cell_peak_ptc_ids cell_peak_ptc_xyz outputFile
% 
% fprintf('Running time of writing results of segmentation: %.2f seconds \n', toc);
%% Cell_backward_forward filtering
% our_dir = 'D:\Point Cloud Processing\Matlab Codes\Bridge decomposition\Version 4\Version 4_';
tic
num_region = length(Cells_Region);
regions_pts = [];
for region_id = 1:num_region
    cell_ptc_ids = vertcat(Cells_Region(region_id).cell.ptc_ids);
    cell_ptc_xyz = OQTR.pts(cell_ptc_ids,1:3);
    cell_ptc_xyz(:,4) = region_id;   
    regions_pts = [regions_pts;cell_ptc_xyz];   
end
% outputFile = strcat(our_dir, strtok(filename,'.'),'_cell_backward_forward_filtering_1.txt');
% outputFile = strcat(our_dir, strtok(filename,'.'),'_cell_forward_region_growing_2.txt');
% outputFile = strcat(our_dir, strtok(filename,'.'),'_cell_boundary_filtering.txt');
% outputFile = strcat(our_dir, strtok(filename,'.'),'_cell_non_plane_filtering.txt');
outputFile = strcat(out_dir, strtok(filename,'.'),'_segments_merging_2_1.txt');
write_txt(outputFile, regions_pts)

clear no_region region_id cell_ptc_ids cell_ptc_xyz regions_pts
fprintf('Running time of writing results: %.2f seconds \n', toc);
% 
% region_id = 11;
% region_ptc_ids = vertcat(Region(region_id).cell.ptc_ids);
% region_ptc_xyz = Tree.pts(region_ptc_ids,:);
% outputFile = strcat(our_dir, 'region_11.txt');
% write_txt(outputFile, region_ptc_xyz)
% clear region_ptc_ids region_ptc_xyz
%%
our_dir = 'D:\Point Cloud Processing\Matlab Codes\Bridge decomposition\Version 4\Version 4_';
tic
num_region = length(cells_region);
regions_pts = [];
for region_id = 1:num_region
    cell_ptc_ids = vertcat(cells_region(region_id).cell.ptc_ids);
    cell_ptc_xyz = OQTR.pts(cell_ptc_ids,1:3);
    cell_ptc_xyz(:,4) = region_id;   
    regions_pts = [regions_pts;cell_ptc_xyz];   
end
% outputFile = strcat(our_dir, strtok(filename,'.'),'_cell_backward_forward_filtering_1.txt');
% outputFile = strcat(our_dir, strtok(filename,'.'),'_cell_forward_region_growing_2.txt');
% outputFile = strcat(our_dir, strtok(filename,'.'),'_cell_boundary_filtering.txt');
% outputFile = strcat(our_dir, strtok(filename,'.'),'_cell_non_plane_filtering.txt');
outputFile = strcat(our_dir, strtok(filename,'.'),'_test_final.txt');
write_txt(outputFile, regions_pts)

clear no_region region_id cell_ptc_ids cell_ptc_xyz regions_pts
fprintf('Running time of writing results: %.2f seconds \n', toc);
%% Write the remain points
tic
flag = inf(size(OQTR.pts,1), 3);
flag(:,1) = 1:size(OQTR.pts,1);
for i = 1:length(cells_plane.cell_ids)
    cell_peak_id = cells_plane.cell_ids(i,:);
    if cell_peak_id(3) == 1
    
        cell_ptc_ids = cells_plane.peak_info(i).ptc_ids;
        cell_peak_id = cells_plane.cell_ids(i,:);
        flag(cell_ptc_ids,2) = cell_peak_id(1);
        flag(cell_ptc_ids,3) = cell_peak_id(2);
    end
end
mask = any(isinf(flag),2);
flag(mask,:) = [];

outputFile = strcat(our_dir, strtok(filename,'.'),'_Remain_Cell_Peak.txt');
write_txt(outputFile, [OQTR.pts(flag(:,1),1:3), flag(:,2:3)])
fprintf('Running time of writing results of point extraction: %.2f seconds \n', toc);
clear cell_peak_id flag mask

%% Write the region of superstructure
our_dir = 'D:\Point Cloud Processing\Matlab Codes\Bridge decomposition\Version 4\Version 4_';

num_component = length(bridge_tree.Component);
for i = 1:num_component
    comp_cell_ptc_ids = vertcat(bridge_tree.Component(i).cell.ptc_ids);
    comp_cell_ptc_xyz = PTC.xyz(comp_cell_ptc_ids,1:3);
    
    outputFile = strcat(our_dir, strtok(filename,'.'),'_Component_', num2str(i), '.txt');
    write_txt(outputFile, comp_cell_ptc_xyz)

end

fprintf('Running time of writing results of point extraction: %.2f seconds \n', toc);

%% Test cell forward region growing
% Region on boundary
mask = ismember(cells.cell_ids(:,1:2), regions_bound_cells_peaks_info(:,1:2), 'rows');
regions_bound_peak_ptc_ids = vertcat(cells.peak_info(mask).ptc_ids);
regions_bound_peak_ptc_xyz = Tree.pts(regions_bound_peak_ptc_ids,:);
outputFile = strcat(our_dir, 'regions_bound_peak_ptc_xyz.txt');
write_txt(outputFile, regions_bound_peak_ptc_xyz)
clear mask regions_bound_peak_ptc_ids regions_bound_peak_ptc_xyz

% Non region
mask = ismember(cells.cell_ids(:,1:2), non_regions_cells_peaks_info(:,1:2), 'rows');
non_regions_peak_ptc_ids = vertcat(cells.peak_info(mask).ptc_ids);
non_regions_peak_ptc_xyz = Tree.pts(non_regions_peak_ptc_ids,:);
outputFile = strcat(our_dir, 'non_regions_peak_ptc_xyz.txt');
write_txt(outputFile, non_regions_peak_ptc_xyz)
clear mask non_regions_peak_ptc_ids non_regions_peak_ptc_xyz

% Points on local surface
outputFile = strcat(our_dir, 'neighbour_cell_peak_ptc_xyz.txt');
write_txt(outputFile, neighbour_cell_peak_ptc_xyz)

outputFile = strcat(our_dir, 'neighbour_check_cell_peak_ptc_xyz.txt');
write_txt(outputFile, neighbour_check_cell_peak_ptc_xyz)

outputFile = strcat(our_dir, 'add_cell_peak_ptc_xyz.txt');
write_txt(outputFile, add_cell_peak_ptc_xyz)

%%



%% Region 2
region_id = 3;
region_ptc_ids = vertcat(cell_region(region_id).cell.ptc_ids);
region_ptc_xyz = Tree.pts(region_ptc_ids,:);
outputFile = strcat(our_dir, 'region_3_r4.txt');
write_txt(outputFile, region_ptc_xyz)
clear region_ptc_ids region_ptc_xyz
%% Road curb
st_cell_ptc_xyz = Tree.pts(st_cell_ptc_ids,:);
outputFile = strcat(our_dir, 'first surface.txt');
write_txt(outputFile, st_cell_ptc_xyz)
clear mask non_regions_peak_ptc_ids non_regions_peak_ptc_xyz

nd_cell_ptc_xyz = Tree.pts(nd_cell_ptc_ids,:);
outputFile = strcat(our_dir, 'second surface.txt');
write_txt(outputFile, nd_cell_ptc_xyz)
clear mask non_regions_peak_ptc_ids non_regions_peak_ptc_xyz



road_curb_cell_ptc_xyz = Tree.pts(road_curb_cell_ptc_ids,:);
outputFile = strcat(our_dir, 'road curb after removing.txt');
write_txt(outputFile, road_curb_cell_ptc_xyz)
clear mask non_regions_peak_ptc_ids non_regions_peak_ptc_xyz


road_curb_cell_ptc_xyz = Tree.pts(road_curb_cell_ptc_ids,:);
outputFile = strcat(our_dir, 'road curb after filtering.txt');
write_txt(outputFile, road_curb_cell_ptc_xyz)
clear mask non_regions_peak_ptc_ids non_regions_peak_ptc_xyz

road_curb_ptc_ids = vertcat(road_curb.cell.ptc_ids);
road_curb_ptc_xyz = Tree.pts(road_curb_ptc_ids,1:3);
outputFile = strcat(our_dir, 'road curb.txt');
write_txt(outputFile, road_curb_ptc_xyz)

% Test the final results in bridge tree

road_curb_ptc_ids = vertcat(bridge_tree. Component(4).cell.ptc_ids);
road_curb_ptc_xyz = Tree.pts(road_curb_ptc_ids,1:3);
outputFile = strcat(our_dir, 'road curb.txt');
write_txt(outputFile, road_curb_ptc_xyz)


%% Parapet extraction

outputFile = strcat(our_dir, 'leaf_cell_cent.txt');
write_txt(outputFile, leaf_cell_cent)


outputFile = strcat(our_dir, 'footpath_polygon.txt');
write_txt(outputFile, footpath_poly)


parapet_cell_ptc_ids = vertcat(Tree.cell_pts(parapet_cell_ids).id);
parapet_cell_ptc_xyz = Tree.pts(parapet_cell_ptc_ids,:);
outputFile = strcat(our_dir, 'footpath_point.txt');
write_txt(outputFile, parapet_cell_ptc_xyz)


outputFile = strcat(our_dir, 'parapet_cell_ptc.txt');
write_txt(outputFile, parapet_cell_ptc_xyz)

outputFile = strcat(our_dir, 'parapet_cell_ptc_dist_filter.txt');
write_txt(outputFile, parapet_cell_ptc_xyz)

outputFile = strcat(our_dir, 'parapet_cell_ptc_footpath_filter.txt');
write_txt(outputFile, parapet_cell_ptc_xyz)



outputFile = 'neighbour_ptc.txt';
write_txt(outputFile, footpath_parapet_cell_ptc_xyz)
% Write results of the parapet

parapet_ptc_ids = vertcat(parapet.cell.ptc_ids);
parapet_ptc_xyz = Tree.pts(parapet_ptc_ids,1:3);
% outputFile = strcat(our_dir, 'parapet.txt');
outputFile = 'parapet.txt';
write_txt(outputFile, parapet_ptc_xyz)

%%
parapet_ptc_xyz = []
for k = 1:length(parapet.cell)
    cell_ptc_ids = parapet.cell(k).ptc_ids;
    cell_ptc_xyz = Tree.pts(cell_ptc_ids,1:3);
    cell_ptc_xyz(:,4) = parapet.cell(k).id(1);
    parapet_ptc_xyz = [parapet_ptc_xyz; cell_ptc_xyz];
    
end
outputFile = 'parapet_cell 1.txt';
write_txt(outputFile, parapet_ptc_xyz)
clear parapet_ptc_xyz  cell_ptc_ids cell_ptc_xyz
%%




level = 1
parapet_ptc_ids = vertcat(final_parapet(level).cell.ptc_ids);
parapet_ptc_xyz = Tree.pts(parapet_ptc_ids,1:3);
outputFile = strcat('parapet_', num2str(level),'.txt');
write_txt(outputFile, parapet_ptc_xyz)


outputFile = 'project_bound_cells.txt';
write_txt(outputFile, parapet_cell_bounds_proj_cent)

outputFile = 'bound_cells.txt';
write_txt(outputFile, parapet_cell_bounds)


mask = ismember(parapet_cell_ids(:,1),cluster_cell_ids );
parapet_ptc_ids = vertcat(parapet.cell(mask).ptc_ids);
parapet_ptc_xyz = Tree.pts(parapet_ptc_ids,1:3);
outputFile = 'test_parapet.txt';
write_txt(outputFile, parapet_ptc_xyz)

%% Cell dist
outputFile = 'cell_idst.txt';
write_txt(outputFile, [parapet_cell_bounds_proj_cent, dist_parapet_road, dist_cell_bin_ids])

%%
level = 2;
cell_ids = vertcat(final_parapet(level).cell.id);
parapet_ptc_xyz = [];
for count = 1:size(cell_ids, 1)
    parapet_ptc_ids = final_parapet(level).cell(count).ptc_ids;
    parapet_ptc_ids(:,2) = cell_ids(count,1);
    cell_ptc_xyz = Tree.pts(parapet_ptc_ids(:,1),1:3);
    parapet_ptc_xyz = [parapet_ptc_xyz; [cell_ptc_xyz,parapet_ptc_ids(:,2)]];
end
    
outputFile = strcat('parapet_cell_ptc', num2str(level),'.txt');
write_txt(outputFile, parapet_ptc_xyz)
clear parapet_ptc_xyz  cell_ptc_xyz parapet_ptc_ids cell_ids


%% Test bottom surface of superstructure
cel_id = [6881,1];
mask = ismember(cells_plane.cell_ids(:,[1,2]),cel_id, 'rows');
cell_ptc_ids = vertcat(cells_plane.peak_info(mask).ptc_ids);
cell_ptc_xyz = Tree.pts(cell_ptc_ids,:);
outputFile = 'road_cell.txt';
write_txt(outputFile, cell_ptc_xyz)
    
cel_id = [6881,2];
mask = ismember(cells_plane.cell_ids(:,[1,2]),cel_id, 'rows');
cell_ptc_ids = vertcat(cells_plane.peak_info(mask).ptc_ids);
cell_ptc_xyz = Tree.pts(cell_ptc_ids,:);
outputFile = 'bottom_cell.txt';
write_txt(outputFile, cell_ptc_xyz)
 

%% Test intermediate
outputFile = 'top.txt';
write_txt(outputFile, Tree.pts(top_surf_cell_ptc_ids,:))
outputFile = 'bottom.txt';
write_txt(outputFile, Tree.pts(bot_surf_cell_ptc_ids,:))
outputFile = 'exclude.txt';
write_txt(outputFile, cell_ptc_xyz)
outputFile = 'intermediate.txt';
write_txt(outputFile, intermediate_cell_ptc_xyz)

outputFile = 'dist_top.txt';
write_txt(outputFile, [cell_ptc_xyz, dist_ptc_top_surface])
outputFile = 'dist_bottom.txt';
write_txt(outputFile, [cell_ptc_xyz, dist_ptc_bot_surface])


%% Write intermediate surface
intermediate_ptc_xyz = [];

for i = 1:length(intermediate_surf.cell)
    cell_ptc_ids = intermediate_surf.cell(i).ptc_ids;
    cell_ptc_xyz = OQTR.pts(cell_ptc_ids,1:3);
    cell_ptc_xyz(:,4) = intermediate_surf.cell(i).ids(1);
    intermediate_ptc_xyz = [intermediate_ptc_xyz;cell_ptc_xyz];
end
% outputFile = strcat(our_dir, strtok(filename,'.'),'_Bridge_component_2.txt');
outputFile = 'intermediate_surface.txt';
write_txt(outputFile, intermediate_ptc_xyz)
%% Write the remaining points
outputFile = 'remaining_ptc.txt';
write_txt(outputFile, OQTR.pts(find(Active_Pts),:))

outputFile = 'dist_cell_traffic.txt';
write_txt(outputFile, [cells_ids(:,1),cells_cent, dist_parapet_road])

%% Segmentation

outputFile = 'segmentation.txt';
write_txt(outputFile, [cell_peak_pts_xyz(ptc_segment_info(:,1),:), ptc_segment_info(:,2)])

outputFile = 'intermediate_surface.txt';
write_txt(outputFile, Tree.pts(intermediate_ptc_ids(flag),:))

%% Check cell
outputFile = 'check_cell.txt';
write_txt(outputFile, Tree.pts(leaf_cell_active_pts_ids,:))

%% Test abutment
% Cell pts
i = 2;
abutment_id = abutment_ids(i);
all_cell_pts_xyz = [];
for k = 1:length(subTree(abutment_id).cell)
    % Retrieve the candidate points of the abutment
    cell_id = subTree(abutment_id).cell(k).ids;
    cell_pts_ids = (subTree(abutment_id).cell(k).ptc_ids);
    cell_pts_xyz = Tree.pts(cell_pts_ids,1:3);
    cell_pts_xyz(:,4) = cell_id(1);
    all_cell_pts_xyz = [all_cell_pts_xyz;cell_pts_xyz];
end

outputFile = 'abutment segment.txt';
write_txt(outputFile, all_cell_pts_xyz)
clear abutment_id cell_id cell_pts_ids cell_pts_xyz all_cell_pts_xyz
%%
leaf_voxel_ids = leafNode(COCT);
all_voxel_ptc_xyz = [];
for k = 1:numel(leaf_voxel_ids)
    
    voxel_id = leaf_voxel_ids(k);

    voxel_ptc_ids = COCT.voxel_ptc_ids(voxel_id).id;
    voxel_ptc_xyz = COCT.pts(voxel_ptc_ids,:);
    voxel_ptc_xyz(:,4) = voxel_id;
    all_voxel_ptc_xyz = [all_voxel_ptc_xyz; voxel_ptc_xyz];
end

outputFile = 'voxel_ptc.txt';
write_txt(outputFile, all_voxel_ptc_xyz)
clear all_voxel_ptc_xyz voxel_ptc_xyz voxel_ptc_ids

%%    region
region_id = 6;
mask = global_regions(:,2) == region_id;
region_voxel_ids = global_regions(mask,1);
all_voxel_ptc_xyz = [];
for k = 1:numel(region_voxel_ids)
    
    voxel_id = region_voxel_ids(k);

    voxel_ptc_ids = COCT.voxel_ptc_ids(voxel_id).id;
    voxel_ptc_xyz = COCT.pts(voxel_ptc_ids,:);
    voxel_ptc_xyz(:,4) = voxel_id;
    all_voxel_ptc_xyz = [all_voxel_ptc_xyz; voxel_ptc_xyz];
end

outputFile = 'region_ptc.txt';
write_txt(outputFile, all_voxel_ptc_xyz)
clear all_voxel_ptc_xyz voxel_ptc_xyz voxel_ptc_ids
%%

all_voxel_ptc_xyz = vertcat(Component_Region(1).voxel.ptc_xyz);
outputFile = 'refine_region_ptc.txt';
write_txt(outputFile, all_voxel_ptc_xyz)

           
%% check voxel
voxel_id = 3814;
voxel_ptc_ids = COCT.voxel_ptc_ids(voxel_id).id;
voxel_ptc_xyz = COCT.pts(voxel_ptc_ids,:);
outputFile = 'a voxel.txt';
write_txt(outputFile, voxel_ptc_xyz)


%%
outputFile = 'interior_voxel.txt';
write_txt(outputFile, interior_voxel_ptc_xyz)


%%

outputFile = 'abutment segment.txt';
write_txt(outputFile, [Tree.pts(ptc_segment_info(:,1),:), ptc_segment_info(:,2)])

%%

outputFile = 'abutment intersection.txt';
write_txt(outputFile, [intersect_line(1:3); intersect_line(1:3) + 10*intersect_line(4:6)])


outputFile = 'abutment outlier.txt';
write_txt(outputFile, [Tree.pts(region_outier_ptc_ids,:)])

%%  Plot voxe vert
all_voxel_vert = [];
for k = 1:numel(neighbour_voxel_ids)
    voxel_id = neighbour_voxel_ids(k);
    voxel_bound = COCT.voxel_bounds(voxel_id,:);
    voxel_vert =  voxel_bound([1 2 3;...
                               4 2 3;...
                               4 5 3;...
                               1 5 3;...
                               1 2 6;...
                               4 2 6;...
                               4 5 6;...
                               1 5 6]);
    voxel_vert(:,4) = voxel_id;
    all_voxel_vert = [all_voxel_vert; voxel_vert];
    
end
outputFile = 'voxel vertices.txt';
write_txt(outputFile, all_voxel_vert)


%% Segmentation for pier

outputFile = 'segmentation.txt';
write_txt(outputFile, [pier_cluster_ptc_xyz(ptc_segment_info(:,1),:), ptc_segment_info(:,2)])


outputFile = 'pier segment.txt';
write_txt(outputFile, [Tree.pts(ptc_segment_info(:,1),:), ptc_segment_info(:,2)])

outputFile = 'voxel.txt';
write_txt(outputFile, interior_voxel_ptc_xyz)

outputFile = '3dbb.txt';
write_txt(outputFile, a)


outputFile = 'new_3dbb.txt';
write_txt(outputFile, bbox_new_verts)

outputFile = 'mid_sect.txt';
write_txt(outputFile, center_sect)


%% Test rotation matrix
outputFile = strcat(our_dir, 'test_rot.txt');
% write_txt(outputFile, rot_region_ptc_xyz)
write_txt(outputFile, (rot_matrix*region_ptc_xyz')')

outputFile = 'target.txt';
write_txt(outputFile, target_ptc_xyz)


outputFile = 'source.txt';
write_txt(outputFile, source_ptc_xyz)


outputFile = 'tersection_line.txt';
write_txt(outputFile, [intersect_line(1:3);intersect_line(1:3)+10*intersect_line(4:6)])


%% Write results of bridge components
all_component_ptc_xyz = [];
level = 5;
for i = 1:length(Bridge(level).Component)
    component_ptc_ids = vertcat(Bridge(level).Component(i).cell.ptc_ids);
    component_ptc_xyz = OQTR.pts(component_ptc_ids,1:3);
    component_ptc_xyz(:,4) = i;
    all_component_ptc_xyz = [all_component_ptc_xyz;component_ptc_xyz];
end
% outputFile = strcat(our_dir, strtok(filename,'.'),'_Bridge_component_2.txt');
outputFile = strcat(our_dir, strtok(filename,'.'),'_Bridge_component_Abutment_2.txt');
outputFile = strcat(our_dir, strtok(filename,'.'),'_Bridge_component_Pier.txt');

write_txt(outputFile, all_component_ptc_xyz)
fprintf('Running time of writing results of point extraction: %.2f seconds \n', toc);
%% Write the points of the substructure
all_component_ptc_xyz = [];
for i = 1:length(bridge_substructure.cell)
    cell_id = bridge_substructure.cell(i).ids;
    component_ptc_ids = bridge_substructure.cell(i).ptc_ids;
    component_ptc_xyz = PTC.xyz(component_ptc_ids,1:3);
    component_ptc_xyz(:,4) = cell_id;
    all_component_ptc_xyz = [all_component_ptc_xyz;component_ptc_xyz];
end
outputFile = strcat(our_dir, strtok(filename,'.'),'_Bridge_Substructure_Candidate_ptc.txt');
write_txt(outputFile, all_component_ptc_xyz)
clear all_component_ptc_xyz cell_id component_ptc_ids component_ptc_xyz
fprintf('Running time of writing results of point extraction: %.2f seconds \n', toc);
%% Write the region of substructure
all_component_ptc_xyz = [];
region_ids = unique(substructure_region_cell_info(:,2));
for i = 1:numel(region_ids)
    mask = substructure_region_cell_info(:,2) == region_ids(i);
    region_cell_ids = substructure_region_cell_info(mask,1);
    for j = 1:numel(region_cell_ids)
        
        mask = ismember(bridge_substructure_cell_ids, region_cell_ids(j));
        region_cell_ptc_ids = bridge_substructure.cell(mask).ptc_ids;
        region_cell_ptc_xyz = PTC.xyz(region_cell_ptc_ids,1:3);
        region_cell_ptc_xyz(:,4) = region_ids(i);
        region_cell_ptc_xyz(:,5) = region_cell_ids(j);
        all_component_ptc_xyz = [all_component_ptc_xyz;region_cell_ptc_xyz];
    end

end
% outputFile = strcat(our_dir, strtok(filename,'.'),'_Bridge_Substructure_region.txt');
outputFile = strcat(our_dir, strtok(filename,'.'),'_Bridge_Substructure_filtering_region.txt');
write_txt(outputFile, all_component_ptc_xyz)
clear all_component_ptc_xyz mask region_cell_ids region_cell_ptc_xyz region_ids
fprintf('Running time of writing results of point extraction: %.2f seconds \n', toc);

%% Compute features for each segments
tic

outputFile = strcat(our_dir, strtok(filename,'.'),'_Bridge_Abutment_1_segmentation.txt');
write_txt(outputFile, [ptc(Segment_Ptc_Info(:,1), 1:3), Segment_Ptc_Info(:,2)])
fprintf('Running time of writing results of point extraction: %.2f seconds \n', toc);    
    
%% Write each surfaces in each ids
tic
all_component_ptc_xyz = [];
component_id = 1;
for i = 1:length(bridge)
    
    sub_component_ids = vertcat(bridge(i).Component.ids);
    for j = 1:numel(sub_component_ids)
        component_ptc_ids = vertcat(bridge(i).Component(j).cell.ptc_ids);
        component_ptc_xyz = PTC.xyz(component_ptc_ids,1:3);
        component_ptc_xyz(:,4) = i;
        component_ptc_xyz(:,5) = component_id;
        all_component_ptc_xyz = [all_component_ptc_xyz;component_ptc_xyz];
        
        % Write for each components
        outputFile = strcat(our_dir, strtok(filename,'.'),'_Bridge_All_Component_', num2str(component_id),'.txt');
        write_txt(outputFile, component_ptc_xyz)
        
        component_id = component_id + 1;

    end
end
outputFile = strcat(our_dir, strtok(filename,'.'),'_Bridge_All_Components.txt');
write_txt(outputFile, all_component_ptc_xyz)

clear all_component_ptc_xyz component_id i j sub_component_ids component_ptc_ids component_ptc_xyz
fprintf('Running time of writing results of point extraction: %.2f seconds \n', toc);

%% Cell region for single region
regions_pts = [];
region_id = 25
region_cell_ids = vertcat(Cell_Region(region_id).cell.id);
for j = 1:size(region_cell_ids,1)       
    cell_ptc_ids = Cell_Region(region_id).cell(j).ptc_ids;
    cell_ptc_xyz = OQTR.pts(cell_ptc_ids,1:3);
    cell_ptc_xyz(:,4) = region_cell_ids(j,1); 
    regions_pts = [regions_pts;cell_ptc_xyz];
end

outputFile = strcat(our_dir, strtok(filename,'.'),'_Forward_region_growing_v2_2.3.txt');
% outputFile = strcat(our_dir, strtok(filename,'.'),'_Boundary_filter_2.txt');
% outputFile = strcat(our_dir, strtok(filename,'.'),'_non_plane_cell_filter.txt');
% outputFile = strcat(our_dir, strtok(filename,'.'),'_region 2_22.txt');
write_txt(outputFile, regions_pts)

clear no_region cell_ptc_ids cell_ptc_xyz regions_pts

%% Cell region for an individual region
num_region = length(Cell_Region);
regions_pts = [];
for region_id = 1:num_region
     
        cell_ptc_ids = vertcat(Cell_Region(region_id).cell.ptc_ids);
        cell_ptc_xyz = OQTR.pts(cell_ptc_ids,1:3);
        cell_ptc_xyz(:,4) = region_id;   
        regions_pts = [regions_pts;cell_ptc_xyz];
end
% outputFile = strcat(our_dir, strtok(filename,'.'),'_backward_forward_filtering_v2_2.txt');
% outputFile = strcat(our_dir, strtok(filename,'.'),'_Forward_region_growing_V2_2.txt');
% outputFile = strcat(our_dir, strtok(filename,'.'),'_Boundary_filtering_V1_1.txt');
outputFile = strcat(our_dir, strtok(filename,'.'),'_Merging_segments_V2_2.txt');
write_txt(outputFile, regions_pts)

clear no_region region_id cell_ptc_ids cell_ptc_xyz regions_pts

% Cell_Region = newRegion
%% Substructure
num_seg = max(substructure_region_cell_info(:,2));
all_region_cell_ptc_xyz = [];
for i = 1: num_seg
    mask = substructure_region_cell_info(:,2) == i;
    region_cell_ids = substructure_region_cell_info(mask,1);
    region_cell_ptc_ids = vertcat(OQTR.cell_pts(region_cell_ids).id);
    region_cell_ptc_xyz = OQTR.pts(region_cell_ptc_ids,1:3);
    region_cell_ptc_xyz(:,4) = i;
    all_region_cell_ptc_xyz = [all_region_cell_ptc_xyz;region_cell_ptc_xyz];
end
    
% outputFile = strcat(our_dir, strtok(filename,'.'),'_sub_segmentation.txt');
outputFile = strcat(our_dir, strtok(filename,'.'),'_filter_sub_segmentation.txt');
write_txt(outputFile, all_region_cell_ptc_xyz)
clear i mask num_seg all_region_cell_ptc_xyz region_cell_ids region_cell_ptc_ids region_cell_ptc_xyz

%% Abutment segmentation
all_region_cell_ptc_xyz = zeros(size(Segment_Ptc_Info,1), 4);
all_region_cell_ptc_xyz(:,1:3) = ptc(Segment_Ptc_Info(:,1),1:3);
all_region_cell_ptc_xyz(:,4) = Segment_Ptc_Info(:,2);

% outputFile = strcat(our_dir, strtok(filename,'.'),'_sub_segmentation.txt');
outputFile = strcat(our_dir, strtok(filename,'.'),'_Abutment_01_segmentation.txt');
write_txt(outputFile, all_region_cell_ptc_xyz)
clear all_region_cell_ptc_xyz


%% pier Column segmentation
all_region_cell_ptc_xyz = zeros(size(Segment_Ptc_Info,1), 4);
all_region_cell_ptc_xyz(:,1:3) = column_candidate_ptc_xyz(Segment_Ptc_Info(:,1),1:3);
all_region_cell_ptc_xyz(:,4) = Segment_Ptc_Info(:,2);

% outputFile = strcat(our_dir, strtok(filename,'.'),'_sub_segmentation.txt');
outputFile = strcat(our_dir, strtok(filename,'.'),'_Pier_01_segmentation_test.txt');
write_txt(outputFile, all_region_cell_ptc_xyz)
clear all_region_cell_ptc_xyz



%% check
all_ptc = [];
for k=1:numel(bound_cell_ids)
    bound_cell_ptc_ids = vertcat(Tree.cell_pts(bound_cell_ids(k)).id);
    bound_cell_ptc_xyz = Tree.pts(bound_cell_ptc_ids, 1:3);
    bound_cell_ptc_xyz(:,4) = bound_cell_ids(k);
    all_ptc = [all_ptc; bound_cell_ptc_xyz];
end
    outputFile = strcat(our_dir, strtok(filename,'.'),'_test_2.txt');
    write_txt(outputFile, all_ptc)
    clear all_ptc

outputFile = strcat(our_dir, strtok(filename,'.'),'_test_colum.txt');
write_txt(outputFile, pier_cap_candidate_ptc_xyz)

outputFile = strcat(our_dir, strtok(filename,'.'),'_test_1.txt');
write_txt(outputFile, segment_ptc_xyz)



outputFile = strcat(our_dir, strtok(filename,'.'),'_test_2.1.txt');
write_txt(outputFile, check_cell_peak_ptc_xyz)

outputFile = strcat(our_dir, strtok(filename,'.'),'_test_2.2.txt');
write_txt(outputFile, current_neighbour_cell_peak_ptc_xyz)

outputFile = strcat(our_dir, strtok(filename,'.'),'_test_2.3.txt');
write_txt(outputFile, add_cell_peak_ptc_xyz)





outputFile = strcat(our_dir, strtok(filename,'.'),'_test_edge_voxels.txt');
write_txt(outputFile, current_seed_edges_voxel_center)

outputFile = strcat(our_dir, strtok(filename,'.'),'_test_edge_voxels.txt');
write_txt(outputFile, [current_check_segment_voxel_center, dist_check_voxel_center_intersection_line])



outputFile = strcat(our_dir, strtok(filename,'.'),'_test_substructure.txt');
write_txt(outputFile, substructure_cell_ptc_xyz)

outputFile = strcat(our_dir, strtok(filename,'.'),'_test_substructure.txt');
write_txt(outputFile, [voxel_features(:,2:4),global_regions(:,2)])


outputFile = strcat(our_dir, strtok(filename,'.'),'_test_segment.txt');
write_txt(outputFile, [COCT.pts(:,1:3),Segment_Ptc_Ids(:,2)])




outputFile = strcat(our_dir, strtok(filename,'.'),'_test_seeding_cell.txt');
write_txt(outputFile, cells_peaks_features(mask,[3,4,5,1,end]))

outputFile = strcat(our_dir, strtok(filename,'.'),'_test_seeding_cell.txt');
write_txt(outputFile, cells_peaks_features(:,[3,4,5,1,end]))

outputFile = strcat(our_dir, strtok(filename,'.'),'_test_22.txt');
write_txt(outputFile, region_current_bound_cell_peak_ptc_xyz)


outputFile = strcat(our_dir, strtok(filename,'.'),'_test_1_1.txt');
write_txt(outputFile, add_cell_peak_ptc_xyz)

outputFile = strcat(our_dir, strtok(filename,'.'),'_test_1_center.txt');
write_txt(outputFile, check_cell_surface_features(1:3))

outputFile = strcat(our_dir, strtok(filename,'.'),'_test_1_dist.txt');
write_txt(outputFile, [current_neighbour_cell_peak_ptc_xyz, dist_neighbour_cell_ptc_cell])




outputFile = strcat(our_dir, strtok(filename,'.'),'_test_1_ptc.txt');
write_txt(outputFile, cell_peak_component_ptc_xyz)



outputFile = strcat(our_dir, strtok(filename,'.'),'_connection_center.txt');
write_txt(outputFile, [vertcat(connect_component_cells.cell.center),dist_cell_bridge_deck, dist_cell_central_line])


outputFile = strcat(our_dir, strtok(filename,'.'),'_connection_center.txt');
write_txt(outputFile, cell_center(flag,:))

outputFile = strcat(our_dir, strtok(filename,'.'),'_connection_center.txt');
write_txt(outputFile, [cell_center,db_dist_cell_bridge_deck])


outputFile = strcat(our_dir, strtok(filename,'.'),'_clustering_result.txt');
write_txt(outputFile, [cell_center,segment_cell(:,2)])




outputFile = strcat(our_dir, strtok(filename,'.'),'_test_1_remain_ptc.txt');
write_txt(outputFile, remain_cell_ptc_xyz)




temp_ptc_ids = vertcat(connect_component_cells.cell.ptc_ids);
temp_ptc = OQTR.pts(temp_ptc_ids,1:3);
outputFile = strcat(our_dir, strtok(filename,'.'),'_connection_ptc.txt');
write_txt(outputFile, temp_ptc)



temp_ptc_ids = connect_component_cells.cell(139).ptc_ids;
temp_ptc = OQTR.pts(temp_ptc_ids,1:3);
outputFile = strcat(our_dir, strtok(filename,'.'),'_connection_ptc.txt');
write_txt(outputFile, temp_ptc)

%%
outputFile = 'test1.txt';
write_txt(outputFile, check_cell_poly_xy)



outputFile = 'test2.txt';
write_txt(outputFile, check_cell_poly_xy_offset)


