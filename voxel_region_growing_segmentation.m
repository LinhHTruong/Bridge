function [ptc_segment_info,Component_Region,COCT] = voxel_region_growing_segmentation(oct_ptc,threshold_voxels, varargin)%,option,out_file)
%% This function is to segment an input data by using region growing based 
% on voxels' feature
%
% INPUT VARIABLES:
%
%   data = [x, y, z]
%   threshold = all thresholds relating to octree generation, segmentation,
%   filtering
%   option = 1: extract segments connectivity
%   option ~= 1: no extract segments connectivity
%
% OUTPUT VARIABLES:
%
%   Segment_Ptc_Ids = [ptc_ids, segment_ids], ptc_ids is local
%   indices
%   Component_Region = data structure 
%       Component_Region.id = segment_i
%       Component_Region.voxel = all voxels in COCT in the region
%       Component_Region.voxel.id = voxel_ids
%       Component_Region.voxel.ptc_ids = ptc_ids in the voxels in the
%       segment
%       Component_Region.voxel.ptc_xyz = ptc coordinate in the voxels in the
%       segment
%
%   COCT = Octree data structure generated from all data points
%   Note: ptc_ids in each voxels were filtered. COCT must used with Component_Region

% STRUCTURE OF FUNCTION
%
%   1. Set up structure for voxels: 
%       - Create octree data structure: empty voxels have the number of the
%       points less than a predefined threshold
%       - Compute features of the voxel where wPCA was used
%   2. Voxel Region growing based residual
%   3. Refinement results of segmentation due to the voxels con occupy
%   points of more than 2 planes
%       - Backward: re-gain the points of the segment (plane) from adjacent
%       segment 
%       - Forward: remove the points in the segment but belong to adjacent
%       segment
%       - The process was applied on the voxels on the boundary of the
%       region, which have adjacent voxels belonging to other regions.
%   4.  Merging non_empty leaf nodel
%       - Points within the empty voxels can be merged to the segments if
%       they located on a local surface of the segment, which was defined
%       as an orthogonal distance to the local surface or voxels
%   5.  Merging the segments
%       - Small segment will be merged to the large segment
%       - Co-planar
%       - Short distance
%       - Connectivity
%
% IMPLEMEMTED FUNCTIONS:
%
%   Octree: Built octree
%   
%
% ASSUMPTION:
%   1. Segment step: 
%       - Voxel having the smallest residual based on points within the 
%       voxel and the fitting plane of these points is considered as the
%       initial voxel (seeding voxels) for the region growing process
%       - Voxels have small deviation in term of normal vectors and
%       residual will be merged into the regions of the seeding voxel
%   2. Voxels represent a local surface as a plane

%
% ABOUTH AUTHOR:
%   Copyright(c) and written by Dr. Linh Truong Hong
%   Affiliation: OLRG, Dept. GRS, CiTG, TU Delft
%   Created day: Jan 03, 2014
%   Last update: May 08, 2019
%
% EXECUTE:
%% Demo
% oct_ptc = pier_cluster_ptc_xyz;block_ptc_xyz;beam_filter_ptc_xyz;
% threshold_voxels = voxel_region_growing_threshold;
% voxel_size = 0.1;
% threshold.angle = 5;
% threshold.residual = 0.01;
% threshold.distance = 0.01;
% option = 1;
%% Constant variables and thresholds
% 
% 
% oct_ptc = comp_pts_xyz;
% threshold_voxels = voxel_region_growing_threshold;
% normal_direction = 'all';
%     time_record = false;

%% 1. Set up structure for voxels
if ~isempty(varargin)
    normal_direction = varargin{1};
    time_record = varargin{2};
else
    normal_direction = 'all';
    time_record = false;
end
if time_record
    tic
end
%%
COCT = OcTree(oct_ptc,'maxSize', threshold_voxels.voxel_size, 'normal_direction', normal_direction);
if time_record
    fprintf('Executing time for octree subdivision: %.2f seconds \n', toc)
end
% Determine a level of octree used for segmentation
if time_record
    tic
end
leaf_voxel_ids = leafNode(COCT);
smallest_voxel_size = (COCT.voxel_bounds(leaf_voxel_ids,[4,5,6]) - COCT.voxel_bounds(leaf_voxel_ids,[1,2,3]));
if min(smallest_voxel_size(:)) <= 0.75*threshold_voxels.voxel_size
    % Extract parent voxels of leaf_nodes
    voxel_depth = max(COCT.depth);
    mask = (COCT.depth == voxel_depth - 1);
    leaf_voxel_ids = find(mask == true);  
    % Remove empty voxel
    leaf_voxel_ids = leaf_voxel_ids(COCT.voxel_prop(leaf_voxel_ids)==1);
end

% Extract leaf node and compute features
Extra_Oct = struct('voxel',[]);
voxel_features = inf(numel(leaf_voxel_ids),8);
voxel_features(:,1) = leaf_voxel_ids;

count_voxel = 1;
for i = 1:numel(leaf_voxel_ids)
    voxel_id = leaf_voxel_ids(i);
    voxel_ptc_ids = COCT.voxel_ptc_ids(voxel_id).id;
    voxel_ptc_xyz = COCT.pts(voxel_ptc_ids,1:3);
    w_interior_voxel_cent = mean(voxel_ptc_xyz, 1);
%     w_ptc_norm = eigenspace(voxel_ptc_xyz, 1);
    [eigVector, eigValue, ~]  = eigenspace(voxel_ptc_xyz,0);
    w_interior_voxel_norm = eigVector(1,:);
        
%     [w_ptc_cent, w_ptc_norm] = wPCA('ptc', voxel_ptc_xyz, 'ini_normal_vector', []);
    dist_ptc_plane = abs(dist_3Dpoints_3Dplane(voxel_ptc_xyz, [w_interior_voxel_cent, w_interior_voxel_norm]));
    residual = sqrt(sum(dist_ptc_plane.^2)/(numel(dist_ptc_plane)));
    if (residual <= threshold_voxels.distance)&&(eigValue(2) >= 0.0001) %Ensure it is not line
        % Filter the outlier points within the voxel
        mask = dist_ptc_plane <= threshold_voxels.distance;
        voxel_inlier_ptc_ids = voxel_ptc_ids(mask);

        if 5 <= numel(voxel_inlier_ptc_ids)
            % Update point within the voxel
            COCT.voxel_ptc_ids(voxel_id).id = voxel_inlier_ptc_ids;

            % Update extra octree
            Extra_Oct.voxel(count_voxel).id = voxel_id;
            Extra_Oct.voxel(count_voxel).ptc_ids = voxel_ptc_ids(~mask);
            count_voxel = count_voxel + 1;

            % Re-compute the features
            voxel_ptc_xyz = voxel_ptc_xyz(mask,1:3);
%             [w_ptc_cent, w_ptc_norm] = wPCA('ptc', voxel_ptc_xyz, 'ini_normal_vector', []);
            w_interior_voxel_cent = mean(voxel_ptc_xyz, 1);
            w_interior_voxel_norm = eigenspace(voxel_ptc_xyz, 1);
            dist_ptc_plane = abs(dist_3Dpoints_3Dplane(voxel_ptc_xyz, [w_interior_voxel_cent, w_interior_voxel_norm]));
            voxel_res = sqrt(sum(dist_ptc_plane.^2)/(numel(dist_ptc_plane)));
            temp_voxel_feature = [w_interior_voxel_cent, w_interior_voxel_norm, voxel_res];
            voxel_features(i,2:end) = temp_voxel_feature;
            clear voxel_ptc_ids voxel_ptc_xyz w_ptc_center w_ptc_normal dist_ptc_plane  voxel_residual temp_voxel_feature%eigen_vector voxel_residual temp_voxel_feature
        end
    else

        Extra_Oct.voxel(count_voxel).id = voxel_id;
        Extra_Oct.voxel(count_voxel).ptc_ids = voxel_ptc_ids;
        count_voxel = count_voxel + 1;
    end
end
% remove any singularity voxel
mask = isinf(voxel_features(:, 4));
voxel_features = voxel_features(~mask,:);
if time_record
    fprintf('Executing time for voxel feature computation: %.2f seconds \n', toc)
end

%% 2. Region growing based on the minimum residual
% Create global region 
if time_record
    tic
end
global_regions = zeros(size(voxel_features,1),2); %[Voxel idx, region idx]
global_regions(:,1) = voxel_features(:,1);
check_voxels = zeros(size(voxel_features,1),2);
check_voxels(:,1) = global_regions(:,1);
check_voxels(:,2) = voxel_features(:,8);
% Region growing process
region_no = 1;
%
while ~isempty(check_voxels)   
    %
    % Find an initial seeding voxel, which has the smallest residual   
    [~,min_ids] = min(check_voxels(:,2));
    % Establish current region
    cur_region = check_voxels(min_ids,1);
    % Establish seeding region
    seed_region = check_voxels(min_ids,1); %Current seed 
    check_voxels = setdiff(check_voxels,check_voxels(min_ids,:),'rows');
    %
    while (~isempty(seed_region))&&(~isempty(check_voxels))
        % Find the first voxel in the searching region
        cur_seed_voxel_id = seed_region(1);
        % Retrive neighbouring voxels of the selected voxel, which are
        % defined in the global reion
        neighbour_voxel_ids = query26_neighbour_voxels(COCT, cur_seed_voxel_id, check_voxels(:,1));
        mask = ismember(neighbour_voxel_ids,global_regions(:,1));
        neighbour_voxel_ids = neighbour_voxel_ids(mask,:);
        %
        if ~isempty(neighbour_voxel_ids)
            % Retrieve features of the seed voxel
            mask = ismember(voxel_features(:,1),cur_seed_voxel_id);
            cur_seed_voxel_cent = voxel_features(mask,2:4);
            cur_seed_voxel_normal = voxel_features(mask,5:7);
            
            % Retrieve features of the neight voxel
            mask = ismember(voxel_features(:,1),neighbour_voxel_ids); 
            neighbour_voxel_cent = voxel_features(mask,2:4);
            neighbour_voxel_normal = voxel_features(mask,5:7);
            
            % Calculate COSIN angle between a seeding voxel and its neighbour voxels
            cosin_neighbour_voxel_normal =  abs(cosine_vectors(cur_seed_voxel_normal, neighbour_voxel_normal));
            dist_neighbour_voxel_cent = abs(dist_3Dpoints_3Dplane(neighbour_voxel_cent, [cur_seed_voxel_cent, cur_seed_voxel_normal]));
            
            % Extracting neighbouring voxels satisfies COSIN angle threshold and distance   
            mask = (cosin_neighbour_voxel_normal >= cos(deg2rad(threshold_voxels.max_angle))&...
                   (dist_neighbour_voxel_cent <= threshold_voxels.distance));
            add_voxel_ids = neighbour_voxel_ids(mask);
            
            % Update a current region by adding voxels satisfied the COSIN angle threshold
            if ~isempty(add_voxel_ids)
                cur_region = union(cur_region,add_voxel_ids ,'rows');
                % Remove added voxels in the examined region out of this region
                mask = ismember(check_voxels(:,1),add_voxel_ids);
                check_voxels = setdiff(check_voxels,check_voxels(mask,:),'rows');
                % Determine added voxel satisfied a RMS residual threshold
                mask = ismember(voxel_features(:,1),add_voxel_ids);
                add_voxel_res = voxel_features(mask, end);   
                mask = add_voxel_res <= threshold_voxels.residual;
                add_voxel_seeding_voxel_ids = add_voxel_ids(mask);
                % Update a seeding region by adding satisfied voxel
                if ~isempty(add_voxel_seeding_voxel_ids)
                    seed_region = union(seed_region,add_voxel_seeding_voxel_ids,'rows');
                end
            end
        end
        % Remove a current seeding point out of the seeding region
        mask = ismember(seed_region,cur_seed_voxel_id);
        seed_region = seed_region(~mask);
    end
    % Update a global region
    if size(cur_region,1) >= 1
        global_regions(ismember(global_regions(:,1),cur_region,'rows'),2) = region_no;  
        region_no = region_no + 1;
    else
        global_regions(ismember(global_regions(:,1),cur_region,'rows'),2) = 0;  
    end
%    
end
if time_record
    fprintf('Executing time for segmentation: %.2f seconds \n', toc)
end
clear mask
clear seeding_region check_voxels current_seeding_voxel_id neighbour_voxel_ids normVectCurVoxel normVectNeighbourVoxels cosinAngle
clear current_region add_voxel_residual region_no min_ids

%% Statistic regions based on the number of voxels in the region and update region order in global_regions
if time_record
    tic
end
num_region = max(global_regions(:,2));
if isempty(num_region)
    ptc_segment_info = [];
    Component_Region = [];
    COCT = [];
    fprintf('SEGMENT IS NOT AVAILABLE\n', toc)
    return;
end


region_count = histcounts(global_regions(:,2), num_region);
region_stat(:, 1) = unique(global_regions(:,2));
region_stat(:, 2) = region_count;

[~, ids] = sort(region_stat(:, 2), "descend");
region_stat = region_stat(ids,:);

% Update region_info
temp_region_ids = inf(size(global_regions, 1),1);
for i = 1:size(region_stat,1)
    mask = ismember(global_regions(:, 2), region_stat(i,1));
    temp_region_ids(mask) = i;
end
%
global_regions(:,2) = temp_region_ids;
clear num_region region_count region_stat ids mask temp_region_ids
if time_record
    fprintf('Executing time for update segment ids: %.2f seconds \n', toc)
end

%% 3. Refinement - Backward + Forward
if time_record
    tic
end
Component_Region = struct('id',[],'voxel',[]);
no_region = max(global_regions(:,2));
% un_region_voxel_ids = setdiff(leaf_voxel_ids,global_regions(:,1));
% region_count = 1;
for i = 1:no_region
    % Retrieve the voxels of the current region
    mask = global_regions(:,2) == i;
    region_voxel_ids = global_regions(mask,1);
    
    % Find voxels on boundary of the region, which has at least 1 neibour outside region
    neighbour_voxels = struct('id',[]);
    flag_bound_voxels = false(numel(region_voxel_ids),1); %[cell_id, 1: bound or 0: interior]
    for j = 1: numel(region_voxel_ids)
        neighbour_voxel_ids = query27VoxelNeighbour(COCT, region_voxel_ids(j),leaf_voxel_ids);
        neighbour_voxels(j).id = neighbour_voxel_ids;
        in_region_voxel_ids = intersect(neighbour_voxel_ids,region_voxel_ids);
        if numel(in_region_voxel_ids) < 8%(1 <= numel(out_region_voxel_ids)) %Bound voxel
            % Although voxel is for 3D, on a plane segment, the interior
            % voxel has 8 neighbour voxels only
            flag_bound_voxels(j) = true;
        end
    end    
    clear j out_region_voxel_ids 

    % Filtering process
    Component_Region(i).id = i;
    count_voxels = 0;
    if numel(region_voxel_ids) == 1
        % For the region only one voxel
        % Retrieve points within the voxel
        cur_voxel_id = region_voxel_ids;
        cur_voxel_ptc_ids = COCT.voxel_ptc_ids(cur_voxel_id).id;
        cur_voxel_ptc_xyz = COCT.pts(cur_voxel_ptc_ids,1:3); 
        
        % Update region 
        count_voxels = count_voxels + 1;
        Component_Region(i).voxel(count_voxels).id = cur_voxel_id;
        Component_Region(i).voxel(count_voxels).ptc_ids = cur_voxel_ptc_ids;
        Component_Region(i).voxel(count_voxels).ptc_xyz = cur_voxel_ptc_xyz;
    else
        % For the region more than one voxels
        for j = 1:numel(region_voxel_ids)

            % Extract points within the voxel belonging to the region
            cur_voxel_id = region_voxel_ids(j);
            cur_voxel_ptc_ids = COCT.voxel_ptc_ids(cur_voxel_id).id;
            
            % Extract points within extract octree (considered as outlier)
            mask = ismember(vertcat(Extra_Oct.voxel.id),cur_voxel_id);
            if any(mask)
                cur_extra_voxel_ptc_ids = Extra_Oct.voxel(mask).ptc_ids;
                cur_voxel_all_ptc_ids = [cur_voxel_ptc_ids;cur_extra_voxel_ptc_ids];
            else
                cur_voxel_all_ptc_ids = cur_voxel_ptc_ids;
            end
            % All ptc_ids in the voxel
            cur_voxel_ptc_xyz = COCT.pts(cur_voxel_all_ptc_ids,1:3);

            % Classify the neighbour voxels
            neighbour_voxel_ids = neighbour_voxels(j).id;
            in_region_voxel_ids = intersect(neighbour_voxel_ids,region_voxel_ids); 
            out_region_voxel_ids = setdiff(neighbour_voxel_ids,region_voxel_ids); 
            in_region_interior_voxel_ids = setdiff(in_region_voxel_ids, region_voxel_ids(flag_bound_voxels));
            
            if any(out_region_voxel_ids) % If the current voxel is the boundary voxels 
                % If no interior voxels -> all neighbour voxels to be used
                if isempty(in_region_interior_voxel_ids)
                    in_region_interior_voxel_ids = in_region_voxel_ids;
                end
                
                % Filter inlier points of the current voxel
                interior_voxel_ptc_ids = vertcat(COCT.voxel_ptc_ids(in_region_interior_voxel_ids).id);
                if numel(interior_voxel_ptc_ids) >= 0.25*numel(cur_voxel_all_ptc_ids)
                    % Estimate the plane through interior voxels
                    interior_voxel_ptc_xyz = COCT.pts(interior_voxel_ptc_ids,1:3);
                    [w_interior_voxel_cent, w_interior_voxel_norm] = wPCA('ptc', interior_voxel_ptc_xyz, 'ini_normal_vector',[]);
                else 
                    [w_interior_voxel_cent, w_interior_voxel_norm] = wPCA('ptc', cur_voxel_ptc_xyz, 'ini_normal_vector',[]);
                end
                
                % Filtering the points of the current voxel
                dist_ptc_plane = dist_3Dpoints_3Dplane(cur_voxel_ptc_xyz, [w_interior_voxel_cent, w_interior_voxel_norm]);
                mask = abs(dist_ptc_plane) <= threshold_voxels.distance;
                cur_voxel_ptc_ids = cur_voxel_all_ptc_ids(mask);
                cur_voxel_ptc_xyz = cur_voxel_ptc_xyz(mask,1:3);
%                 else
%                     w_interior_voxel_cent = [];
%                     w_interior_voxel_norm = [];
%                     cur_voxel_ptc_ids = [];
%                 end

                % Backward: Update current voxels for the current region
                if ~isempty(cur_voxel_ptc_ids)
                    count_voxels = count_voxels + 1;                    
                    % Update region
                    Component_Region(i).voxel(count_voxels).id = cur_voxel_id;
                    Component_Region(i).voxel(count_voxels).ptc_ids = cur_voxel_ptc_ids;
                    Component_Region(i).voxel(count_voxels).ptc_xyz = cur_voxel_ptc_xyz;
                end

                % Forward: Add voxels in the out_region: voxels may be in un_region or region
                for k = 1:numel(out_region_voxel_ids)
                    % Retrieve points in the voxel considering that inlier
                    out_region_voxel_id = out_region_voxel_ids(k);
                    out_region_voxel_ptc_ids = COCT.voxel_ptc_ids(out_region_voxel_id).id;
                    
                    % Retrieve points in the voxel considering that outlier
                    mask = ismember(vertcat(Extra_Oct.voxel.id),out_region_voxel_id);
                    if any(mask)
                        out_region_extra_voxel_ptc_ids = Extra_Oct.voxel(mask).ptc_ids;
                        out_region_all_voxel_ptc_ids = [out_region_voxel_ptc_ids; out_region_extra_voxel_ptc_ids];
                    else
                        out_region_all_voxel_ptc_ids = out_region_voxel_ptc_ids;
                    end
                    % All vptc_ids in the voxel
                    out_region_voxel_ptc_xyz = COCT.pts(out_region_all_voxel_ptc_ids,1:3);
                    
                    % Extract inlier points to be merge to the large region
                    dist_ptc_plane = dist_3Dpoints_3Dplane(out_region_voxel_ptc_xyz, [w_interior_voxel_cent, w_interior_voxel_norm]);
                    mask = abs(dist_ptc_plane) <= threshold_voxels.distance;
                    out_region_voxel_inlier_ptc_ids = out_region_all_voxel_ptc_ids(mask);
                    out_region_voxel_inlier_ptc_xyz = out_region_voxel_ptc_xyz(mask,1:3);
%                     out_region_voxel_inlier_ptc_ids = out_region_voxel_ptc_ids(mask);
%                     out_region_voxel_inlier_ptc_xyz = out_region_voxel_ptc_xyz(mask,1:3);
                    
                    % Extract outlier points to be kept
                    out_region_voxel_outlier_ptc_ids = out_region_all_voxel_ptc_ids(~mask);
                    
                    % Update new voxel associated points for the current region
                    if ~isempty(out_region_voxel_inlier_ptc_ids)
                        count_voxels = count_voxels + 1;
                        Component_Region(i).voxel(count_voxels).id = out_region_voxel_ids(k);
                        Component_Region(i).voxel(count_voxels).ptc_ids = out_region_voxel_inlier_ptc_ids;
                        Component_Region(i).voxel(count_voxels).ptc_xyz = out_region_voxel_inlier_ptc_xyz;
                    end

                    % Update for the orginal octree                        
                    COCT.voxel_ptc_ids(out_region_voxel_id).id = out_region_voxel_outlier_ptc_ids;
                end

            else
                % Remain: Update current voxels for the current region
                % Filter the voxels itself
                mask = voxel_features(:,1) == cur_voxel_id;
                cur_voxel_cent = voxel_features(mask,2:4);
                cur_voxel_normal = voxel_features(mask,5:7);
                dist_ptc_plane = dist_3Dpoints_3Dplane(cur_voxel_ptc_xyz, [cur_voxel_cent, cur_voxel_normal]);
                mask = abs(dist_ptc_plane) <= threshold_voxels.distance;
                
                % Update structure
                count_voxels = count_voxels + 1;
                Component_Region(i).voxel(count_voxels).id = cur_voxel_id;
                Component_Region(i).voxel(count_voxels).ptc_ids = cur_voxel_all_ptc_ids(mask);
                Component_Region(i).voxel(count_voxels).ptc_xyz = cur_voxel_ptc_xyz(mask,1:3);
            end
        end
    end
end
if time_record
    fprintf('Running time of the backward and forward filtering: %.2f seconds \n', toc);
end

%% Remove empty region
if time_record
    tic
end
no_region = length(Component_Region);
temp_Component_Region = struct('id',[],'voxel',[]);
region_count = 0;
for i=1:no_region
    if ~isempty(Component_Region(i).voxel)
        region_count = region_count + 1;
        temp_Component_Region(region_count).id = region_count;
        temp_Component_Region(region_count).voxel = Component_Region(i).voxel;
    end
end
%
clear Component_Region
Component_Region = temp_Component_Region;
clear temp_Component_Region
if time_record
    fprintf('Running time for removing empty voxel: %.2f second\n',toc)
end
%% 5. Merging small regions into the larger region
% Merging the voxels in the small region into the large region based on 
% the distance between the points in the voxel in the small region to the
% fitting surface of the points in the voxels belong to the large region
% Compute the region features
if time_record
    tic
end
num_regions = length(Component_Region);
comp_region_features = inf(num_regions,8);%[id, cent_x, cent_y, cent_z, nx, ny, nz, no.ptc]
comp_region_features(:,1) = vertcat(Component_Region.id);
for i = 1:num_regions
    region_ptc_xyz = vertcat(Component_Region(i).voxel.ptc_xyz);
    if size(region_ptc_xyz,1) >= 5
        region_cent = mean(region_ptc_xyz, 1);
        region_normal = eigenspace(region_ptc_xyz, 1);
        comp_region_features(i,2:end) = [region_cent, region_normal, size(region_ptc_xyz,1)];
    end
end
mask = isinf(comp_region_features(:,5));
comp_region_features = comp_region_features(~mask,:);
[~, sort_ids] = sort(comp_region_features(:,end),'descend');
comp_region_features = comp_region_features(sort_ids,:);
clear mask sort_ids num_regions region_ptc_xyz region_cent region_normal 

%% Merging process
% Determine the searching window size
if min(smallest_voxel_size(:)) <= threshold_voxels.voxel_size
    window_size = 5*min(smallest_voxel_size(:)); % The real voxel size = 2*min(smallest_voxel_size(:))
else
    window_size = 2.5*min(smallest_voxel_size(:));
end
% Prealloocation
check_region_ids = comp_region_features(:,1);
update_Region = zeros(size(comp_region_features,1),2);
update_Region(:,1) = comp_region_features(:,1);

% 
count = 1;
while ~isempty(check_region_ids)
%     numel(check_region_ids)
    if size(check_region_ids,1) == 1
        mask = ismember(update_Region(:,1), check_region_ids, 'rows');
        update_Region(mask, 2) = count;
        break
    end
    
    % Update seeding region
    seed_region_ids = check_region_ids(1);
    
    % Update region 
    mask = ismember(update_Region(:,1), seed_region_ids, 'rows');
    update_Region(mask, 2) = count;
    
    % Update check_patch_ind
    check_region_ids(1) = [];
    
    % Searching merging region
    while ~isempty(seed_region_ids)
        
        % Extract features of the seed_region
        cur_seed_region_id = seed_region_ids(1);
        cur_seed_region_voxel_ids = vertcat(Component_Region(cur_seed_region_id).voxel.id);
        mask = ismember(comp_region_features(:,1), cur_seed_region_id);
        cur_seed_region_features = comp_region_features(mask,:);
%         cur_seed_region_cent = comp_region_features(mask,2:4);
%         cur_seed_region_normal = comp_region_features(mask,5:7);
%         
                   
        % Extract target region features
        tg_region_ids = check_region_ids;
        mask = ismember(comp_region_features(:,1), tg_region_ids);
        tg_region_features = comp_region_features(mask,:);
%         tg_region_cent = comp_region_features(mask,2:4);
%         tg_region_normal = comp_region_features(mask,5:7);

        % Check similarity: direction -> normal deviation; 
        cosin_seed_tg_regions = cosine_vectors(cur_seed_region_features(5:7), tg_region_features(:,5:7));
        mask = abs(cosin_seed_tg_regions) >= cos(deg2rad(threshold_voxels.max_angle));
        tg_region_ids = tg_region_ids(mask,:);

        % Check local conunity
        if ~isempty(tg_region_ids)
            
            
            % Check for each part of region
            flag_tg_regions = false(numel(tg_region_ids),1);    
            for k = 1:numel(tg_region_ids)
                % Retrieve the voxels within the 
                tg_region_id = tg_region_ids(k);
                tg_region_voxel_ids = vertcat(Component_Region(tg_region_id).voxel.id);
        
                % searching neighbour voxels in a target region for current seed region
                seed_tg_voxel_ids = windowQueryVoxelNeighbour(COCT, cur_seed_region_voxel_ids, tg_region_voxel_ids, window_size);
                tg_seed_voxel_ids = windowQueryVoxelNeighbour(COCT, tg_region_voxel_ids, cur_seed_region_voxel_ids, window_size);
                
                % Retrieve points within neighbour region
                if ~isempty(seed_tg_voxel_ids)&&~isempty(tg_seed_voxel_ids)
                    % Points in seed region
                    mask = ismember(cur_seed_region_voxel_ids,tg_seed_voxel_ids); 
                    seed_voxel_ptc_xyz = vertcat(Component_Region(cur_seed_region_id).voxel(mask).ptc_xyz);
                    
                    % Points in target region
                    mask = ismember(tg_region_voxel_ids,seed_tg_voxel_ids);
                    tg_voxel_ptc_xyz = vertcat(Component_Region(tg_region_id).voxel(mask).ptc_xyz);
                    
                    % Estimate the surface
                    surface_ptc_xyz = union(seed_voxel_ptc_xyz, tg_voxel_ptc_xyz, 'rows');
                    [surface_cent, surface_norm] = wPCA('ptc', surface_ptc_xyz, 'ini_normal_vector',[]);
                
                     % Compute inlier ratio for points in seed region
                    dist_seed_voxel_ptc_xyz = dist_3Dpoints_3Dplane(seed_voxel_ptc_xyz, [surface_cent, surface_norm]);
                    mask = abs(dist_seed_voxel_ptc_xyz) <= threshold_voxels.distance;
                    ratio_seed_voxel_inlier_ptc_xyz = sum(mask)>= 0.75*size(seed_voxel_ptc_xyz,1);
                    
                    % Compute inlier ratio for points in target region
                    dist_tg_voxel_ptc_xyz = dist_3Dpoints_3Dplane(tg_voxel_ptc_xyz, [surface_cent, surface_norm]);
                    mask = abs(dist_tg_voxel_ptc_xyz) <= threshold_voxels.distance;
                    ratio_tg_voxel_inlier_ptc_xyz = sum(mask)>= 0.75*size(tg_voxel_ptc_xyz,1);
                    
                    if ratio_seed_voxel_inlier_ptc_xyz&&ratio_tg_voxel_inlier_ptc_xyz
                        flag_tg_regions(k) = true;
                    end
                    
                end
            end
            
            % Update merging region
            merge_region_ids = tg_region_ids(flag_tg_regions);
            if ~isempty(merge_region_ids)
                % Merging - Update new region id
                mask = ismember(update_Region(:,1), merge_region_ids, 'rows');
                update_Region(mask, 2) = count;
                % Update seeding region
                mask = ismember(comp_region_features(:,1),merge_region_ids);
                add_region_num_ptc = comp_region_features(mask,end);
                mask = add_region_num_ptc >= 25;
                add_region_ids = merge_region_ids(mask);
                if ~isempty(add_region_ids)
                    seed_region_ids = union(seed_region_ids, add_region_ids);
                end
                % Remove target region to be merge
                mask = ismember(check_region_ids,merge_region_ids);
                check_region_ids(mask) = [];
%                 check_region_ids = setdiff(check_region_ids,merge_region_ids);
            end
        end
        % remove the seed region to be check
        seed_region_ids = setdiff(seed_region_ids, cur_seed_region_id);        
    end  
    count = count + 1;
end
if time_record
    fprintf('Running time for merging the small segments: %.2f seconds \n', toc);
end
clear all_region_cell_ids
%% Update 
if time_record
    tic
end
newComponent_Region = struct('id',[],'voxel',[]);
for i = 1:length(unique(update_Region(:,2)))
    mask = update_Region(:,2) == i;
    % Retrive information from current re
    old_region_ids = update_Region(mask,1);
    % Update a new region
    newComponent_Region(i).id = i;
    newComponent_Region(i).voxel = Component_Region(old_region_ids(1)).voxel;
    if numel(old_region_ids) > 1
        count_id = length(newComponent_Region(i).voxel);   
        for j = 2:numel(old_region_ids)
            num_region_voxels = length(Component_Region(old_region_ids(j)).voxel);
            for k = 1:num_region_voxels
                count_id = count_id + 1;
                newComponent_Region(i).voxel(count_id).id = Component_Region(old_region_ids(j)).voxel(k).id;
                newComponent_Region(i).voxel(count_id).ptc_ids = unique(Component_Region(old_region_ids(j)).voxel(k).ptc_ids);
                newComponent_Region(i).voxel(count_id).ptc_xyz = Component_Region(old_region_ids(j)).voxel(k).ptc_xyz;
            end
        end
    end 
end
%
clear Component_Region
Component_Region = newComponent_Region;
clear newComponent_Region
if time_record
    fprintf('Running time for updating a new segment: %.2f seconds \n', toc);
end
%% Export data
num_regions = length(Component_Region);
ptc_segment_info = [];
for i = 1:num_regions
    region_ptc_ids = unique(vertcat(Component_Region(i).voxel.ptc_ids));
    region_ptc_ids(:,2) = i;
    ptc_segment_info = [ptc_segment_info;region_ptc_ids];
end

% Write output file
%     
% if (option == 1)|(strcmp(option, 'y'))|(strcmp(option, 'yes'))    
%     write_txt(out_file, [COCT.pts(:,1:3),Segment_Ptc_Ids(:,2)])
% end


