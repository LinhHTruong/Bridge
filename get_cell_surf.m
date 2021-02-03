function [cell_ptc_ids, cell_surf] = get_cell_surf(struct_comp, cell_id, option)
%% This function is to retrieve points' indices and surfaces of the cell
% Input:
% Output:
% Demo:
%%
% Get the top surface directly
cell_ptc_ids = vertcat(struct_comp(cell_id).ptc_ids);
cell_surf = vertcat(struct_comp(cell_id).surface_features);
% Select the upper surface
if size(cell_surf,1) > 1
    if strcmp(option, 'top')
        [~, mask] = min(cell_surf(:,3));
    else
        [~, mask] = min(cell_surf(:,3));
    end
    cell_surf = cell_surf(mask,:);
end

end