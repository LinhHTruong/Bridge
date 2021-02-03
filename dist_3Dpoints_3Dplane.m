function d = dist_3Dpoints_3Dplane(ptc, plane)
%DISTANCEPOINTPLANE Signed distance betwen 3D point and plane
%
%   D = distancePointPlane(POINT, PLANE)
%   Returns the euclidean distance between point POINT and the plane PLANE,
%   given by: 
%   ptc : [x0 y0 z0]
%   PLANE : [x0 y0 z0 nx, ny, nz]
%   D     : scalar  
%   
%   See also:
%   planes3d, points3d, intersectLinePlane
%


%   HISTORY


% normalized plane normal
plane_n = normalizeVector3d(plane(:,4:6));
plane_center = plane(:,1:3);
% Uses Hessian form, ie : N.p = d
% I this case, d can be found as : -N.p0, when N is normalized
d = -sum(bsxfun(@times, plane_n, bsxfun(@minus, plane_center, ptc)), 2);