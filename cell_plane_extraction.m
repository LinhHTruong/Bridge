function [Cells_Plane, Cells_nonPlane] = cell_plane_extraction(Tree, cell_ids, threshold)
%% Thi fuction is to extract the cells contain the point cloud of the planes
%(either vertical or horizontal)
% The cell will be splitted into a planar point and non-planar points

% Input:
% Output:
% Developed by: Dr. Linh Truong-Hong, ORSG, Dept. GRS, TUDelft
%               
% %% Demo
% Tree = OQTR;
% cell_ids = leaf_cell_ids;
% threshold = THRESHOLD;
% plane_norm = 'oz'
% peak_loc = 'top' 


%% Check input variable
% Create initial data structures
tic
Cells_Plane = struct('cell_ids',[],'peak_info',[]);
Cells_nonPlane = struct('cell_ids',[],'ptc_ids', []);

% Calculate bandwidth
bandwidth = 0.5*threshold.min_dist_2_planes;

% Establish thresholds

count_plane_cell = 0;
for i=1:numel(cell_ids)
    % Extract points in a cell
    cell_id = cell_ids(i);
    cell_ptc_ids = Tree.cell_pts(cell_id).id;
    cell_ptc_xyz = Tree.pts(cell_ptc_ids,1:3);
    
    % Extract the points in the cell after removing coincide points in the
    % vertical direction
    [~,ia,~] = unique(cell_ptc_xyz(:,1:2), 'rows');
    cell_filter_ptc_xyz = cell_ptc_xyz(ia,:);
    
   % Kernel density estimation
    cell_filter_ptc_range = [min(cell_filter_ptc_xyz(:,3)), max(cell_filter_ptc_xyz(:,3))];
    no_kde_pt = max(100,ceil(diff(cell_filter_ptc_range)/bandwidth));
    [fi,zi,~] = ksdensity(cell_filter_ptc_xyz(:,3),'npoints',no_kde_pt,'bandwidth',bandwidth,'Kernel','epanechnikov');
%     mask = fi < 0.05;
%     i
%     fi = fi(mask)
    peak_shape = peak_shape_width(zi, fi);
    mask = peak_shape(:,3) + peak_shape(:,2) > 0;
    peak_shape = peak_shape(mask,:);
    % Sort from heigh to low
    [~, mask] = sort(peak_shape(:,1), 'descend');
    peak_shape = peak_shape(mask,:);
    clear cell_filter_ptc_xyz
    peak_shape(:,4) = peak_shape(:,2) + peak_shape(:,3);
    peak_shape(:,[2,3]) = [];
     
%     close all
%     hold all
%     plot(zi,fi,'r-')

    % Compute feature of the points within the peak shape (bell)
    peak_count = 1;  
    for j = 1:size(peak_shape,1)
        
        % Extract the points within the shape
        mask = (peak_shape(j,1) - peak_shape(j,2)/2 <= cell_ptc_xyz(:,3))&(cell_ptc_xyz(:,3) <= peak_shape(j,1) + peak_shape(j,2)/2);
        plane_ptc_ids = cell_ptc_ids(mask);
        plane_ptc_xyz = cell_ptc_xyz(mask,:);
        
        % Update Cell Plane 
        if numel(plane_ptc_ids) >= threshold.min_num_pts
            % Filter the surface
%             [w_ptc_center, w_ptc_normal] = wPCA('ptc', plane_ptc_xyz, 'angle_threshold', 0.5*pi/180, 'ini_normal_vector', []);
%             dist_ptc_plane = dist_3Dpoints_3Dplane(plane_ptc_xyz, [w_ptc_center, w_ptc_normal]);
%             mask  = abs(dist_ptc_plane) <= threshold.max_distance;
%             plane_ptc_ids = plane_ptc_ids(mask);
%             plane_ptc_xyz = plane_ptc_xyz(mask, :);
%             if numel(plane_ptc_ids) >= threshold.min_num_pts
                % Compute a surface features
                plane_cent = mean(plane_ptc_xyz,1);
                [plane_eigen_vects, ~, plane_res] = eigenspace(plane_ptc_xyz, 4);
                plane_normal = plane_eigen_vects(1,:);

                % Update the data structure
                count_plane_cell = count_plane_cell + 1;
                Cells_Plane.cell_ids(count_plane_cell,:) = [cell_id, peak_count, 1]; %status = 1: active
                Cells_Plane.peak_info(count_plane_cell).ptc_ids = plane_ptc_ids;
                Cells_Plane.peak_info(count_plane_cell).peaks_features = [plane_cent, plane_normal, plane_res];
                peak_count = peak_count + 1;
                clear plane_eigen_vects plane_cent plane_normal plane_res
%             end
        end
        clear mask peak_ptc_ids peak_ptc_xyz
    end    
    
end
 
% Update non_planar_cell
non_plane_cell_ids = setdiff(cell_ids, unique(vertcat(Cells_Plane.cell_ids(:,1))));
for i = 1:numel(non_plane_cell_ids)
    Cells_nonPlane.cell_ids(i,1) = non_plane_cell_ids(i);
    cell_ptc_ids = Tree.cell_pts(non_plane_cell_ids(i)).id;
    Cells_nonPlane.ptc_ids(i,1).id = cell_ptc_ids; 
    clear cell_ptc_ids
end
clear non_plane_cell_ids 
   
fprintf('Running time for extracting points on a plane: %.2f seconds \n', toc);
