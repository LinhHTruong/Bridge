function proj_ptc = proj_point_on_3Dline(ptc, line)
    %%This function is to project 3D points onto 3D line
    % Input:
    %    ptc           : [Nx3] - x, y, z coordinates
    %    line          : [1x6] - P0 (x, y, z), t (tx, ty, tz)
    % Output:

    %    proj_ptc      : [Nx3]

    % Demo:
    %    ptc = [0, 3, 0]
    %    line = [0,0,0, 1,1,0]


    %%
    p0 = line(1:3);
    vl = line(4:end);
    w = ptc - p0;
    b = sum(bsxfun(@times, w, vl),2)/dot(vl, vl);
    proj_ptc = p0 + b.*vl;
   