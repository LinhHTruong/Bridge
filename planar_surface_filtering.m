function ptc_seg_info = planar_surface_filtering(segs_ptc, ptc_seg_info, check_option, struct_threshold, threshold, time_record)
%% This function is to filter a small planar surfaces derived from the segmentation
% This is done through the minimum dimensions of the structure
% Input:
%       segs_ptc                : [Nx3] - x, y, z coordinates
%       ptc_seg_info            : [Nx2] - [pts_ids, seg_ids]
%       check_option            : Various options:  Full for all dimensions;
%                                                   Width for the cross-section width; 
%                                                   Length for the cross-section length
% Output:
%       ptc_seg_info            : [Nx2] - [pts_ids, seg_ids]
% Demo:
% segs_ptc = cell_peak_pts_xyz;
% ptc_seg_info = ptc_segment_info;
% check_option = 'length';
%% Calculate the features of each segment
if time_record
    tic
end
seg_ids = unique(ptc_seg_info(:,2));
seg_features = inf(numel(seg_ids), 3); %[seg_id, b, h]
for i = 1:numel(seg_ids)
    % Retrieve the points within the segment
    seg_id = seg_ids(i);
    mask = ismember(ptc_seg_info(:,2), seg_id);
    seg_pts_ids = ptc_seg_info(mask,1);
    seg_pts_xyz = segs_ptc(seg_pts_ids,1:3);
    
    % Compute the feature of each segment
    seg_plane_normal = eigenspace(seg_pts_xyz,1);
        
    % Compute dimensions
    rot_euler_angle = vrrotvec(seg_plane_normal,threshold.nz);
    rot_matrix = vrrotvec2mat(rot_euler_angle);
    rot_seg_pts_xyz = (rot_matrix*seg_pts_xyz')';
    seg_mbb = min2DBoundingBox(rot_seg_pts_xyz(:,1:2)');
    clear rot_euler_angle rot_matrix rot_surf_ptc_xyz 

    % Assign surface features
    seg_features(i,:) = [seg_id, seg_mbb.short_edge, seg_mbb.long_edge];
end

%% Filtering
if strcmp(check_option, 'full')
    mask = (max(struct_threshold.section_height, struct_threshold.section_height) <= seg_features(:,2))&(struct_threshold.section_length <= seg_features(:,3));
elseif strcmp(check_option, 'width')
    mask = (max(struct_threshold.section_height, struct_threshold.section_height) <= seg_features(:,2));
elseif strcmp(check_option, 'length')
    mask = (0.05 <= seg_features(:,2))&(struct_threshold.section_length <= seg_features(:,3)); % to avoid the line
else
    error('Check option by %s does not support', check_option);
    return;
end
seg_features = seg_features(mask,:);

%% Update the segment
if ~isempty(seg_features) 
    % Update ptc_segment_info
    mask = ismember(ptc_seg_info(:,2), seg_features(:,1));
    ptc_seg_info(~mask, :) = [];
else
    ptc_seg_info = [];
end

if time_record
    fprintf('Running time for filtering the planar surface based primary beam geometry: %.2f seconds\n',toc);
end
