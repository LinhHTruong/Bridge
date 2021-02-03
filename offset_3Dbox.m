function bbox_new_verts = offset_3Dbox(box_verts, offset)
% This function is to offset the edges of the polygon, which is based on
% the vector between two edges
% Input:
%       poly_verts      : [Nx2]
%       offset          : scalar
% Output:
%       poly_new_verts  : [Nx2]
% Demo:
% box_verts = a
% offset = struct_threshold.section_width


% Set up link-vertices
linked_vert = [1, 2, 4, 5;...
               2, 1, 3, 6;...
               3, 2, 4, 7;...
               4, 1, 3, 8;...
               5, 1, 6, 8;...
               6, 5, 7, 2;...
               7, 6, 8, 3;
               8, 5, 7, 4];

% Compute vectors of edges connected to the vertices
bbox_new_verts = zeros(size(box_verts, 1),3);
for i = 1: size(box_verts, 1)
    
    % Calculate edge vectors
    edge_vects = box_verts(linked_vert(i,1),:) - box_verts(linked_vert(i,2:end),:) ;
    edge_vects = normalizeVector3d(edge_vects);
    
    % The directional vector
    offset_direction = mean(edge_vects,1);
    
    % Determine the new vertices
    bbox_new_verts(i,:) = box_verts(linked_vert(i,1),:) + offset*offset_direction;

end