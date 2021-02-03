function region_features = cal_region_features(Tree, region)
%% This function is to calculate the features of the region, which is preassumed as the planar
% Input:
%       Tree                    : Data structure stored entire data

% Output:

% Developed by: Linh Truong-Hong, Dept. GRS, TU Delft

% Demo:
% region = ptc_segment_info;

% tic

% Preallocation Features
if isstruct(region)
    region_ids = [1:length(region)]';
else
    region_ids = unique(region(:,2));
end
region_features = struct('ids', [], 'center',[],'normal',[],'tangent',[], 'residual',[],'bounding_box',[]);
region_features.ids = region_ids;


% Compute the features of the region
nz = [0., 0., 1.0];
for i = 1:numel(region_ids)
    % Retrieve points within the region
    region_id = region_ids(i);
    if isstruct(region)
        region_ptc_ids = vertcat(region(region_id).cell.ptc_ids);
    else
        mask = region(:,2) == region_id;
        region_ptc_ids = region(mask,1);
    end
    region_ptc_xyz = Tree.pts(region_ptc_ids, 1:3);
    region_features.center(i,:) = mean(region_ptc_xyz, 1);
    region_eigenspace = eigenspace(region_ptc_xyz, 0);
    region_features.normal(i,:) = region_eigenspace(1,:);
    region_features.tangent(i,:) = region_eigenspace(3,:);
    % Compute smooth of the surface/region
    dist_ptc_region = dist_3Dpoints_3Dplane(region_ptc_xyz, [region_features.center(i,:), region_features.normal(i,:)]);
    mean_dist = mean(abs(dist_ptc_region));
    residual_dist = sqrt(sum(dist_ptc_region.^2)/numel(dist_ptc_region));
    region_features.residual(i,:) = [mean_dist,residual_dist];
    
    % Estimate the bounding box
    rot_euler_angle = vrrotvec(region_eigenspace(1,:),nz);
    rot_matrix = vrrotvec2mat(rot_euler_angle);
    rot_region_ptc_xyz = (rot_matrix*region_ptc_xyz')';
    rot_region_mbb = min2DBoundingBox(rot_region_ptc_xyz(:,1:2)');
    region_features.bounding_box(i,:) = [rot_region_mbb.area, rot_region_mbb.long_edge,rot_region_mbb.short_edge];
   
end
% fprintf('Running time of computing features of the regions: %.2f seconds \n', toc);
clear num_regions region_ptc_xyz rregion_eigenspace dist_ptc_region mean_dist residual_dist min_boudingbox