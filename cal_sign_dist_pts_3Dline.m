function d = cal_sign_dist_pts_3Dline(point, line)
%DISTANCEPOINTLINE3D Euclidean distance between 3D point and line
%
%   d = distancePointLine3d(POINT, LINE);
%   Returns the distance between point POINT and the line LINE, given as:
%   POINT : [x0 y0 z0]
%   LINE  : [x0 y0 z0 dx dy dz]
%   D     : (positive) scalar  

% cf. Mathworld (distance point line 3d)  for formula
d = bsxfun(@rdivide, vectorNorm3d( ...
        vectorCross3d(line(:,4:6), bsxfun(@minus, line(:,1:3), point)) ), ...
        vectorNorm3d(line(:,4:6)));

    point = [2,1,1];
    line = [1,1,1,1,1,1]
d = bsxfun(@rdivide,sum(vectorCross3d(line(:,4:6), bsxfun(@minus, line(:,1:3), point)),2),vectorNorm3d(line(:,4:6)))

vectorNorm3d