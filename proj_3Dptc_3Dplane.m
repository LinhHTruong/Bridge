function proj_ptc = proj_3Dptc_3Dplane(ptc, plane)
%PROJPOINTONPLANE Return the orthogonal projection of a point on a plane
%
%   proj_ptc = proj_3Dptc_3Dplane(ptc, plane)
%   ptc = [Nx3]
%   plane = [1x6] : [cent_x, cent_y, cent_z, nx, ny, nz)

% Unpack the planes into origins and normals
surface_center = plane(:,1:3);
surface_normal_vector = plane(:,4:6);
surface_normal_vector = surface_normal_vector./norm(surface_normal_vector);

c = dot(surface_normal_vector, surface_center);
a = sum(bsxfun(@times,ptc,surface_normal_vector), 2) - c;
proj_ptc = ptc - bsxfun(@times, a, surface_normal_vector);
    
    
    
    
    
    