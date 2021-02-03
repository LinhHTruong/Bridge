function ptc_sign = ptc_sign_distance(ptc, line)
%% This function is to determine the sign of the points regarding to the line
% Input:
%     ptc         : Nx3
%     line        : 1x6
% Output:
%     ptc_sign    : Nx1 (-1, 1)
% Demo:
%     ptc
%     line
%%    
% Extract information
line_center = line(1:3);
line_tangent = line(4:6);

% Convert to xp plane
line_center(3) = 0.0;
line_tangent(3) = 0.0;
line_tangent = line_tangent./norm(line_tangent);

% Project points to the line
ptc_sign = ones(size(ptc,1),1);
proj_ptc = proj_point_on_3Dline(ptc, [line_center, line_tangent]);
proj_ptc_ptc = ptc - proj_ptc;
ptc_sign(2:end) = sign(sum(bsxfun(@times,proj_ptc_ptc(1,:),proj_ptc_ptc(2:end,:)),2));%proj_ptc_ptc(1,:)- normal
    