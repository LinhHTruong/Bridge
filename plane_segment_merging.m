function segs_info = plane_segment_merging(segs_ptc, segs_info, struct_threshold, threshold, time_record)

% This function is used to merge the segment from the points or voxel-based
% region growing method. This is only applied for a linear structural
% components likely columns or beams
% Input:
%       segs_ptc                : [Nx3]
%       segs_info               : [Nx2] - [pts_ids, seg_ids]
% Output:

% Development:
%     Linh Truong-Hong
% Demo:

%   segs_ptc = cell_peak_pts_xyz;%block_ptc_xyz;
%   segs_info = ptc_segment_info;
%   comp_name = 'beam'
%   debug = 'false'
if time_record
    tic
end

%% 1. Compute the features of the region
% nz = [0.,0.,1.0];
seg_ids = unique(segs_info(:,2));
segs_features = inf(numel(seg_ids),8);%[id, cent_x, cent_y, cent_z, nx, ny, nz, no.ptc]
segs_features(:,1) = seg_ids;
for i = 1:numel(seg_ids)    
    % Retrieve the points within the region
    mask = segs_info(:,2) == seg_ids(i);
    seg_ptc_ids = segs_info(mask,1);
    seg_ptc = segs_ptc(seg_ptc_ids,1:3);
    
    % Compute segment features
    if size(seg_ptc,1) >= 20
        seg_cent = mean(seg_ptc, 1);
        seg_normal = eigenspace(seg_ptc, 1);
        segs_features(i,2:end) = [seg_cent, seg_normal, size(seg_ptc,1)];
    end
end

mask = isinf(segs_features(:,5));
segs_features = segs_features(~mask,:);
[~, sort_ids] = sort(segs_features(:,end),'descend');
segs_features = segs_features(sort_ids,:);
clear i mask seg_ids seg_ptc_ids seg_ptc sort_ids

%% Merging process

% Prealloocation
check_seg_ids = segs_features(:,1);
update_seg_ids = zeros(size(segs_features,1),2);
update_seg_ids(:,1) = segs_features(:,1);

% 
count = 1;
while ~isempty(check_seg_ids)
%     numel(check_region_ids)
    if size(check_seg_ids,1) == 1
        mask = ismember(update_seg_ids(:,1), check_seg_ids, 'rows');
        update_seg_ids(mask, 2) = count;
        break
    end
    
    % Select initial region
    ini_seed_segment_id = check_seg_ids(1);
    % Update region 
    mask = ismember(update_seg_ids(:,1), ini_seed_segment_id, 'rows');
    update_seg_ids(mask, 2) = count;
    
    % Update seeding region
    seed_seg_ids = ini_seed_segment_id;
    
    % Update check_patch_ind
    check_seg_ids(1) = [];
    
    %% Searching merging region
    while ~isempty(seed_seg_ids)
        
        % Extract features of the ref_region
        seed_region_id = seed_seg_ids(1);
        mask = ismember(segs_features(:,1), seed_region_id);
        seed_seg_features = segs_features(mask,:);
                
        % Extract sample region features
        mask = ismember(segs_features(:,1), check_seg_ids);
        tg_segs_features = segs_features(mask,:);

        % Check similarity: direction -> normal deviation
        cosin_segs = cosine_vectors(seed_seg_features(5:7), tg_segs_features(:,5:7));
        mask = abs(cosin_segs) >= cos(deg2rad(2.0*threshold.max_angle));
        tg_segs_features = tg_segs_features(mask,:);

        % Check distance
        if ~isempty(tg_segs_features)
            dist_segs = dist_3Dpoints_3Dplane(tg_segs_features(:,2:4), seed_seg_features(2:7));
            mask = abs(dist_segs) <= 2.0*threshold.max_distance;
            merge_region_features = tg_segs_features(mask,:);
            
            % Overlap
            flag_merging_region = false(numel(merge_region_features(:,1)),1);
            
            % Retrieve points within the source segment
            mask = ismember(segs_info(:,2),seed_region_id);
            seed_region_ptc_ids = segs_info(mask, 1);
            seed_region_ptc_xyz = segs_ptc(seed_region_ptc_ids,1:3);
            for i = 1:size(merge_region_features,1)
                
                % Retrieve points within the target region
                mask = ismember(segs_info(:,2),merge_region_features(i,1));
                tg_region_ptc_ids = segs_info(mask, 1);
                tg_region_ptc_xyz = segs_ptc(tg_region_ptc_ids,1:3);
                
                % Searching overlap
                % Searching neigbour for source segment
                [sr_ptc_ids, tg_ptc_ids] = cloud_cloud_overlap(seed_region_ptc_xyz, tg_region_ptc_xyz, struct_threshold);
%                 
%                 [neighbour_ptc_ids, neighbour_dist] = knnsearch(tg_region_ptc_xyz, seed_region_ptc_xyz, 'k', 1000);
%                 mask = neighbour_dist <= struct_threshold.section_width;
%                 tg_ptc_ids = unique(neighbour_ptc_ids(mask)); % neighbour of seed_ptc in tg_region
%                 mask = all(~mask,2);
%                 sr_ptc_ids = find(mask == 0);
                
                if (numel(sr_ptc_ids) > 20)&&(numel(tg_ptc_ids) > 20)
                    % overlap surface
                    sr_ptc_xyz = seed_region_ptc_xyz(sr_ptc_ids,1:3);
                    tg_ptc_xyz = tg_region_ptc_xyz(tg_ptc_ids,1:3);
                                       
                    % Check directions of the sub-datasets: 
                        overlap_ptc_xyz = union(sr_ptc_xyz, tg_ptc_xyz, 'rows'); 
                        overlap_center = mean(overlap_ptc_xyz, 1);
                        overlap_normal = eigenspace(overlap_ptc_xyz, 1);

                        % Check inlier
                        dist_sr_ptc_overlap_surface = dist_3Dpoints_3Dplane(sr_ptc_xyz, [overlap_center,overlap_normal]);
                        mask = abs(dist_sr_ptc_overlap_surface) <= 2.0*threshold.max_distance;
                        sr_inlier_ptc_no = sum(mask) >= 0.5*numel(sr_ptc_ids);

                        dist_tg_ptc_overlap_surface = dist_3Dpoints_3Dplane(tg_ptc_xyz, [overlap_center,overlap_normal]);
                        mask = abs(dist_tg_ptc_overlap_surface) <= 2.0*threshold.max_distance;
                        tg_inlier_ptc_no = sum(mask) >= 0.5*numel(tg_ptc_ids);

                        if sr_inlier_ptc_no&&tg_inlier_ptc_no
                            flag_merging_region(i) = true;
                        end
%                     end
                end
            end
            
            % Update merging region
            merge_region_ids = merge_region_features(flag_merging_region,1);
        
            % Update merging 
            if ~isempty(merge_region_ids)
                % Merging - Update new region id
                mask = ismember(update_seg_ids(:,1), merge_region_ids, 'rows');
                update_seg_ids(mask, 2) = count;
                % Update seeding region
                seed_seg_ids = union(seed_seg_ids, merge_region_ids);
                % Remove target region to be merge
%                 check_seg_ids = setdiff(check_seg_ids, merge_region_ids);
                mask = ismember(check_seg_ids, merge_region_ids);
                check_seg_ids(mask) = [];
            end
        end

        % remove the seed region to be check
        seed_seg_ids = setdiff(seed_seg_ids, seed_region_id);
    end  
    count = count + 1;
end

%% Update new segments

temp_seg_ids = inf(size(segs_info,1),1);
segs_new_ids = unique(update_seg_ids(:,2));
for i = 1:numel(segs_new_ids)
    
    % Retrive information from current re
    mask = update_seg_ids(:,2) == segs_new_ids(i);
    seg_old_id = update_seg_ids(mask,1);
    seg_new_id = unique(update_seg_ids(mask,2));
    
    % Update a new region
    mask = ismember(segs_info(:,2), seg_old_id);
    temp_seg_ids(mask) = seg_new_id; 
end

% Replace the old ids by the new ids
segs_info(:,2) = temp_seg_ids;

% Remove sigularity segment
mask = isinf(segs_info(:,2));
if any(mask)
    segs_info(mask,:) = [];
end
if time_record
    fprintf('Running time for updating a new segment: %.2f seconds \n', toc);
end

end % End of the function

%% Find an overlap
function [seed_neighbour_ptc_ids, target_neighbour_ptc_ids] = cloud_cloud_overlap(seed_pts_xyz, target_pts_xyz, struct_threshold)
%% This function is to find points' indicies in an overlap area defined by the minimum bounding box (mBB) plus an offset
% Input:
% Output:
% Demo:
% seed_pts_xyz;
% target_pts_xyz;


%% Rotation the data
nz = [0., 0., 1.0];
% Estimate the surface of the seed points
seed_surf_cent = mean(seed_pts_xyz,1);
seed_surf_normal = eigenspace(seed_pts_xyz,1);

% Project seed_pts and target_pts on the seed surface
% Seed data
seed_proj_pts_xyz = proj_3Dptc_3Dplane(seed_pts_xyz, [seed_surf_cent, seed_surf_normal]);
% Tagrget data
target_proj_pts_xyz = proj_3Dptc_3Dplane(target_pts_xyz, [seed_surf_cent, seed_surf_normal]);

% Project onto xy plane
rot_euler_angle = vrrotvec(seed_surf_normal, nz);
rot_matrix = vrrotvec2mat(rot_euler_angle);
seed_rot_ptc_xyz = (rot_matrix*seed_proj_pts_xyz')';
target_rot_ptc_xyz = (rot_matrix*target_proj_pts_xyz')';


%% Determine the bounding box
% Minimum bounding box
seed_mbb = min2DBoundingBox(seed_rot_ptc_xyz(:,1:2)');
target_mbb = min2DBoundingBox(target_rot_ptc_xyz(:,1:2)');

% Get veetices of the polygon
seed_poly = seed_mbb.vertices;
target_poly = target_mbb.vertices;

% Create offset polygon
seed_poly_offset = polygon_offset(seed_poly, struct_threshold.section_width);
target_poly_offset = polygon_offset(target_poly, struct_threshold.section_width);

%% Find points in the polygon
% For seed
[in, on] = inpolygon(seed_rot_ptc_xyz(:,1), seed_rot_ptc_xyz(:,2), target_poly_offset(:,1), target_poly_offset(:,2));
seed_neighbour_ptc_ids = find(in | on);
% For target
[in, on] = inpolygon(target_rot_ptc_xyz(:,1), target_rot_ptc_xyz(:,2), seed_poly_offset(:,1), seed_poly_offset(:,2));
target_neighbour_ptc_ids = find(in | on);

end %End of the function


