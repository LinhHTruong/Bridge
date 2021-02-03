function ptc_local_id = road_curb_ransac_filter(pts, threshold, debug)
%% This function is to extract planes of structures (column, beam) by using ransace
% Input:
%   ptc                 : [Nx3] - x, y, z coordinates

% Output:
%   ptc_seg_info        : [Nx2] - [ptc_ids, segment_ids]    
% Demo:
%   pts = road_curb_cell_ptc_xyz;
%
if nargin <= 2
    
   debug = false;
end
    
if debug
    tic
end

%% Initial point cloud
pts = pointCloud(pts);

%% Preallocation
Model = struct('id',[], 'model', [], 'ptc_ids',[]);
flag = true;
model_no = 1;
input_ptc = pts;
input_ptc_ids = 1:pts.Count;
%
while flag
    % Estimate the model
    [model,inlier_ptc_local_ids,outlier_ptc_local_ids] = pcfitplane(input_ptc, 0.5*threshold.max_distance);
    
    % Assess the model
    inlier_ptc_ids = input_ptc_ids(inlier_ptc_local_ids);
    cosin_normal = cosine_vectors(model.Normal, threshold.nz);
    if (abs(cosin_normal) <= cos(deg2rad(2.0*threshold.max_angle)))&&(numel(inlier_ptc_ids) > threshold.min_num_pts)
        % Update data structure
        Model(model_no).id = model_no;
        Model(model_no).model = model;
        Model(model_no).ptc_ids = inlier_ptc_ids;
        model_no = model_no + 1;
    end
    
    % Set up data for a next iteration
    if numel(outlier_ptc_local_ids) >= threshold.min_num_pts
        input_ptc_ids = input_ptc_ids(outlier_ptc_local_ids);
        input_ptc = select(input_ptc,outlier_ptc_local_ids);
    else
        flag = false;
    end
    clear large_region_id
end

%% Gathering points in segment
ptc_seg_info = zeros(pts.Count,2);
ptc_seg_info(:,1) = 1:pts.Count;
for i = 1:length(Model)
    % Retrieve points of the model
    ptc_ids = Model(i).ptc_ids;
    mask = ismember(ptc_seg_info(:,1), ptc_ids);
    ptc_seg_info(mask,2) = i;
end

% Extract the largest segment
mask = isinf(ptc_seg_info(:, 2)) | (ptc_seg_info(:,end) == 0);
ptc_seg_info(mask,:) = [];

if ~isempty(ptc_seg_info)
    num_region = max(ptc_seg_info(:, 2));
    region_count = histcounts(ptc_seg_info(:, 2), num_region);
    region_stat(:, 1) = unique(ptc_seg_info(:,end));
    region_stat(:, 2) = region_count;

    [~, max_id] = max(region_stat(:, 2));
    large_region_id = region_stat(max_id,1);
    mask = ptc_seg_info(:, 2) == large_region_id;
    ptc_local_id = ptc_seg_info(mask,1);
else
    ptc_local_id = [];
end

if debug
    fprintf('Running time of the filtering based ransac: %.2f seconds \n', toc);
end   