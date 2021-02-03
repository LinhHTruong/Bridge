function d = dist_3Dpoints_3Dline(ptc, line, varargin)
%% Compute distances from the points to lines
% Input:
%     ptc         - [Nx3]
%                             
%     line        - [1x6] - cent_x, cent_y, cent_z, tx, ty, tz
%     varargin    - [1x3] the vector perpendicular to the line to get sign distance 
% Output:
%     d           - [Nx1]
% Demo: 
%% Compute the distance (possitive)
proj_ptc = proj_point_on_3Dline(ptc, line);
proj_ptc_ptc = ptc - proj_ptc;
d = sqrt(sum(bsxfun(@power, proj_ptc_ptc, 2),2));

if numel(varargin) == 1
    tt = varargin{:};
    if numel(tt) == 3
        ptc_sign = sign(sum(bsxfun(@times,tt,proj_ptc_ptc),2));%proj_ptc_ptc(1,:)- normal
        d = d.*ptc_sign;
    else
        warning('Input must be a 3D vector');
        return;
    end
elseif numel(varargin) > 1
    warning('Number of input variables excceeded');
    return;

end
