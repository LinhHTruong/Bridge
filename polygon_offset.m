function poly_new_verts = polygon_offset(poly_verts, offset)
% This function is to offset the edges of the polygon, which is based on
% the vector between two edges
% Input:
%       poly_verts      : [Nx2]
%       offset          : scalar
% Output:
%       poly_new_verts  : [Nx2]
% Demo:
% poly_verts = seed_region_poly
% offset = struct_threshold.section_width


% Check polygon
if any(poly_verts(1,:) == poly_verts(end,:))
    % A close polygon
    poly_verts = poly_verts(1:end-1,:);
end

% Add head and tail vertices
poly_verts = [poly_verts(end,:);poly_verts;poly_verts(1,:)];

% Compute a pair of vectors at each vertice
poly_vert_edge_vect_1 = poly_verts(2:end-1,:) - poly_verts(1:end-2,:);
poly_vert_edge_vect_2 = poly_verts(2:end-1,:) - poly_verts(3:end,:);

% Offset direction at each vertice
offset_direction = (normalizeVector3d(poly_vert_edge_vect_1) + normalizeVector3d(poly_vert_edge_vect_2))/2;
%     offset_direction = normalizeVector3d(offset_direction);
% Determine the new vertices
poly_new_verts = poly_verts(2:end-1,:) + offset*offset_direction;


poly_new_verts(end+1,:) = poly_new_verts(1,:);
end