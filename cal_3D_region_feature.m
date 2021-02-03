function region_features = cal_3D_region_feature(Tree, region_info)
%% This function is to filter the segments of the column
% Input:
% Output:
% Demo:

% region_info = ptc_segment_info;
%% Compute the features of the column
if isstruct(region_info)
    region_ids = [1:length(region_info)]';
else
    region_ids = unique(region_info(:,2));
end
region_features = struct('ids', [], 'center',[], 'normal',[],'tangent',[], 'edge_length',[]);

% Calculate the cluster height
cluster_ptc_xyz = Tree.pts(region_info(:,1), 1:3);
cluster_height = max(cluster_ptc_xyz(:,3)) - min(cluster_ptc_xyz(:,3));

% Calculate the feature
count = 1;
for i = 1:numel(region_ids)
    % Retrieve points within the region
    region_id = region_ids(i);
    if isstruct(region_info)
        region_ptc_ids = vertcat(region_info(region_id).cell.ptc_ids);
    else
        mask = region_info(:,2) == region_id;
        region_ptc_ids = region_info(mask,1);
    end
    region_ptc_xyz = Tree.pts(region_ptc_ids, 1:3);
    
    if numel(region_ptc_ids) < 20
        continue
    end
    
    % Compute principal direction
    region_eigenspace = eigenspace(region_ptc_xyz, 0);

    % Estimate 3D bounding box
%     [~,bb_corner,~,~,bb_edge_vectors, bb_edge_length] = min3Dboundingbox(region_ptc_xyz(:,1),region_ptc_xyz(:,2),region_ptc_xyz(:,3),'v',3);
    [~,~,~,~,~, bb_edge_length] = min3Dboundingbox(region_ptc_xyz(:,1),region_ptc_xyz(:,2),region_ptc_xyz(:,3),'v',3);
    
    % Assign the features
    if bb_edge_length(end) < 0.5*cluster_height
        continue;
    end
    region_features.ids(count,1) = region_id;
    region_features.center(count,:) = mean(region_ptc_xyz, 1);
    region_features.normal(count,:) = region_eigenspace(1,:);
    region_features.tangent(count,:) = region_eigenspace(3,:);
    region_features.edge_length(count,:) = bb_edge_length(:);
%     region_features.verts(count).id = bb_corner;
    count = count + 1;
end

%%
% 
% 
% 
% %%
% edge_list = [1, 2; 2, 3; 3, 4; 4, 1;...
%              5, 6; 6, 7; 7, 8; 8, 5;...
%              1, 5; 2, 6; 3, 7; 4, 8];
% surface_list = [1, 2, 3, 4;...
%                 5, 6, 7, 8;...
%                 1, 2, 6, 5;...
%                 2, 3, 7, 6;...
%                 3, 4, 8, 7;...
%                 4, 1, 5, 8];
% mid_edge = 0.5*(bbox_new_verts(edge_list(:,1),:) + bbox_new_verts(edge_list(:,2),:))
% mid_surface = 0.25*(bbox_new_verts(surface_list(:,1),:) + bbox_new_verts(surface_list(:,2),:) + ...
%                     bbox_new_verts(surface_list(:,3),:) + bbox_new_verts(surface_list(:,4),:))
% 
% 
% outputFile = 'mid_sect.txt';
% write_txt(outputFile, [mid_edge; mid_surface])
%                 
