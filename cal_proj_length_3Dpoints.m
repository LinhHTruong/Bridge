function [end_proj_pts, proj_length]  = cal_proj_length_3Dpoints(ptc, line, threshold_distance)
% This function is to find the end points and its length when projecting
% The points within the certain distance

%     Input:
%     Output
%     Demo:

if ~isempty(threshold_distance)
    % Compute distance from the vertices to the line
    dist_ptc_line = dist_3Dpoints_3Dline(ptc, line);

    % Find closest points
    mask = dist_ptc_line <= threshold_distance;
    ptc = ptc(mask,:);
end

% Find the project length
if size(ptc,1) >= 2
    proj_ptc = proj_point_on_3Dline(ptc, line);
    [~, ~, end_proj_pts, ~, ~] = get_end_points(proj_ptc);
    end_proj_pts = reshape(end_proj_pts, [3,2])';
    proj_length = norm(diff(end_proj_pts,1,1));

else
    end_proj_pts = [];
    proj_length = 0;
end
       
end % end of the function