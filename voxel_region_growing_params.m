function voxel_region_growing_threshold = voxel_region_growing_params(section_type, edge_length, threshold_angle, sampling_step)
    % This function is to compute the parameters for the voxel_base region
    % growing
    % Demo:
    % section_type = 1 - planar
	% section_type = 0 - Curve
    % edge_length = col_surface_info(9);
    % threshold_angle = threshold.max_angle;
    % sampling_step = threshold.sampling_step;
            % Voxel-based region growing

    % Preallocation
    voxel_region_growing_threshold = struct('voxel_size',[],'residual',[], 'distance',[], 'max_angle',[], 'min_num_pts', []);
    
    % Compute the parameters 
    if section_type
        % Rectangle section
        % Voxel size
        voxel_size = max(0.05, edge_length/2);
        % Residual
        res = sampling_step;
        % Distance threshold
        dist = 2.0*sampling_step;
        % Angle threshold
        max_angle = threshold_angle;
    else
        scale = [0.75, 1.5];
        voxel_size = sqrt(2)*edge_length*sin(deg2rad(threshold_angle)/2);
        if voxel_size < 10*sampling_step
            % Re-adjust the voxel size
            voxel_size = 10*sampling_step;
            threshold_angle = rad2deg(2*asin(voxel_size/(sqrt(2)*edge_length)));
        end
        
        % Residual
        res = scale(2)*edge_length*(1 - cos(deg2rad(threshold_angle)/2));
        
        % Distance threshold
        dist = scale(2)*2*edge_length*sin(deg2rad(threshold_angle))*sin(deg2rad(threshold_angle)/2);

        % Angle threshold
        max_angle = scale(2)*threshold_angle;
    end
    % Assign the parameters
    voxel_region_growing_threshold.voxel_size = voxel_size;
    voxel_region_growing_threshold.residual = res;
    voxel_region_growing_threshold.distance = dist;
    voxel_region_growing_threshold.max_angle = max_angle;
    % The minimum number of the points within the voxel
    voxel_region_growing_threshold.min_num_pts = 5;
    
end