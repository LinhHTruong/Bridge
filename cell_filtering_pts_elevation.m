function decision = cell_filtering_pts_elevation(pts_z, struct_threshold)
%% This function is to filter the cell may/may not containing the points of the substructure
% If the cell has a segment with its length larger than the min_length of the
% structural component, the cell represent a substructure
%
%   Input:
%   Output:
%   Demo:
% pts_z = leaf_cell_active_pts_z;

%% Compute the gap between two consecutive points
pts_z = sort(pts_z);
diff_ptc_z = diff(pts_z);

%% Determine the gap
% Preallocation of a decision
decision = false;
% Check the cells
mask = 0.5*struct_threshold.column_min_height <= diff_ptc_z;
if any(mask)
    % Get loc_ids at the end of the continuous segment (start of the gap)
    end_gap_loc_ids = find(mask);
    % Get loc_ids at the start of the continuous segment (end of the gap)
    start_gap_loc_ids = end_gap_loc_ids + 1; 

    % Adding head and tails of each loc_ids
    start_gap_loc_ids = union(1,start_gap_loc_ids);
    end_gap_loc_ids = union(end_gap_loc_ids, numel(pts_z));

    % Adjust loc_ids in the range of the data
    start_gap_loc_ids = start_gap_loc_ids(start_gap_loc_ids <= numel(pts_z));
    
    % Remove coincide loc_ids (both start and end loc_ids have the same id)
    mask = start_gap_loc_ids ~= end_gap_loc_ids;
    start_gap_loc_ids = start_gap_loc_ids(mask);
    end_gap_loc_ids = end_gap_loc_ids(mask);

    % Calculate a length of the continuous segment
    seg_length = pts_z(end_gap_loc_ids) - pts_z(start_gap_loc_ids);
    
    % Get a decision
    if max(seg_length) >= struct_threshold.column_min_height
        decision = true;
%     else
%         decision = false;
    end
else
    if abs(pts_z(end) - pts_z(1)) >= struct_threshold.column_min_height
        decision = true;
    end
end

end %End of the function