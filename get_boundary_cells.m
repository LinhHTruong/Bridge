function [region_bound_cell_ids, region_interior_cell_ids] = get_boundary_cells(TREE, region_cell_ids)
%% This function is to find the cells on the boundary of the 2D regions.
% Input:
%       TREE    : Data structure of the cell
%       region_cell_ids : a list of cell ids in the region
% Output
% Demo
%   TREE = OQTR
%   region_cell_ids = region_cell_ids


% Find cells on boundary of the region, which has less than 8 neibour
% region
% neighbour_cell = struct('id',[]);
cell_bound_info = zeros(numel(region_cell_ids),2); %[cell_id, 1: bound or 0: interior]
cell_bound_info(:,1) = region_cell_ids;

for j = 1: numel(region_cell_ids)
    neighbour_cell_ids = Query_Neighbour_Cells(TREE, region_cell_ids(j),region_cell_ids);
%     neighbour_cell(j).id = neighbour_cell_ids;
    if (numel(neighbour_cell_ids) >=1 ) && (numel(neighbour_cell_ids) < 8) %Bound cells
        cell_bound_info(j,2) = 1; %Bound cells
    end
end

region_bound_cell_ids = cell_bound_info(logical(cell_bound_info(:,2)),1);
region_interior_cell_ids = cell_bound_info(~logical(cell_bound_info(:,2)),1);
