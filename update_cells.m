function cells_plane = update_cells(Tree, cells_plane, cells_region, un_assigned_segs_ids,... 
                        super_structure, superstructure_comp_name,threshold, time_record)
% function cells_plane = update_cells(Tree, cells_plane, cells_region, super_structure, superstructure_comp_name, un_assigned_segs_ids, threshold, time_record)

%% This function is to update the status of cells in the cell_planes. 

% Cells assinged as the supstructure are deactivated. 
% Cells above defined surfaces of the superstructure are also deactivated.
% Input:
%     Tree                      - 2D quadtree
%     cells_plane               - Planes in each 2D cells
%     cells_region              - Regions of cells
%     super_structure           - Data structure stored components of the bridge
%     superstructure_comp_name  - A level of the surface components stored. 1: Top surface of a superstructure; 2: Bottom surface; 3: Immediate surface
%     un_assigned_segs_ids      - Segments were extracted from a cell-based region growing but have not assigned as the surface component at the certain level
%     threshold                 - A filtering threshold
    
% Output:
%     cells_plane             -  Planes in each 2D cells
%% Demo:
% Tree = OQTR;
% cells_plane = Cells_Plane;
% cells_region = Cells_Region;
% super_structure = SuperStructure
% superstructure_comp_name = {};
% un_assigned_segs_ids = un_connect_segs_ids;
% threshold = THRESHOLD;
%% Recover all cells in Regions but not recognize that the bridge components
% This is becaue
if time_record
    tic
end
for i = 1:numel(un_assigned_segs_ids)
    % Retrieve the cells in the segment
    un_assigned_seg_id = un_assigned_segs_ids(i);
    seg_cell_ids = vertcat(cells_region(un_assigned_seg_id).cell.id);
    
    % Recover the cells and their features as original (After cell plane extraction)
    for j = 1:size(seg_cell_ids,1)
        cur_seg_cell_peak_id = seg_cell_ids(j,:);
        mask = ismember(cells_plane.cell_ids(:,[1,2]),cur_seg_cell_peak_id, 'rows'); 
        if any(mask)
            % Update status for the cell_peak
            cells_plane.cell_ids(mask,end) = 1;
            cell_peak_ptc_ids = cells_plane.peak_info(mask).ptc_ids;
            
            % Retrieve ptc_ids within the cell-peak in the region
            cur_seg_cell_peak_ptc_ids = cells_region(un_assigned_seg_id).cell(j).ptc_ids;
            if ~isequal(cell_peak_ptc_ids,cur_seg_cell_peak_ptc_ids) 
                % Compute features of the cell-peak
                cell_peak_ptc_ids = union(cell_peak_ptc_ids,cur_seg_cell_peak_ptc_ids);
                cell_peak_ptc_xyz = Tree.pts(cell_peak_ptc_ids,1:3);

                cell_peak_cent = mean(cell_peak_ptc_xyz,1);
                [cell_peak_eigen_vectors, ~, cell_peak_surface_res] = eigenspace(cell_peak_ptc_xyz, 4);
                cell_peak_normal = cell_peak_eigen_vectors(1,:);

                % Update the features of the cell-peak
                cells_plane.peak_info(mask).ptc_ids = cell_peak_ptc_ids;
                cells_plane.peak_info(mask).peaks_features = [cell_peak_cent, cell_peak_normal, cell_peak_surface_res];
                clear cell_peak_ptc_ids cell_peak_ptc_xyz cell_peak_center cell_peak_eigen_vectors cell_peak_normal cell_peak_surface_residual
            end
        end
    end
end


%% Remove all sigularity of the cell: points within the peak was merged to the region assigned as the bridge components
new_cells_plane = struct('cell_ids',[],'peak_info',[]);
count = 1;
for i = 1:length(cells_plane.cell_ids)
    cell_peak_features = cells_plane.peak_info(i).peaks_features;
    if ~isinf(cell_peak_features)
        new_cells_plane.cell_ids(count,:) = cells_plane.cell_ids(i,:);
        new_cells_plane.peak_info(count).ptc_ids = cells_plane.peak_info(i).ptc_ids;
        new_cells_plane.peak_info(count).peaks_features = cells_plane.peak_info(i).peaks_features;
        count = count + 1;
    end
end
cells_plane = new_cells_plane;
clear i new_Cells cell_peak_features count

%% Update the status of the cells: The cells above the road surface, footpath are deactived 
% Retrieve cell_peak features from Cells
cells_peaks_ids = cells_plane.cell_ids;
inactive_cell_peak_ids = cells_peaks_ids(cells_peaks_ids(:,3) == 0,:);
cell_peaks_features = vertcat(cells_plane.peak_info.peaks_features);

% Turn status of all cell peak to active
cells_plane.cell_ids(:,end) = 1;
for i = 1:length(super_structure.Component)
    % Retrieve the component names
    comp_name = super_structure.Description{i};
    if any(contains(comp_name,superstructure_comp_name, 'IgnoreCase', true))
        continue;
    end
    % Retrieve cell_ids in the components
    comp_cell_ids = vertcat(super_structure.Component(i).cell.ids);

    for j = 1:size(comp_cell_ids,1)
        % Sewarch cell peaks
        comp_cell_id = comp_cell_ids(j,:);

        if comp_cell_id(2) ~= 0
            % The cells assigned as the road surface and footpaths
            mask = (cells_plane.cell_ids(:,1) == comp_cell_id(1))&(cells_plane.cell_ids(:,2) <= comp_cell_id(2));
        else
            % For the cells as parapets
%             cell_ptc_ids = super_structure.Component(i).cell(j).ptc_ids;
%             cell_cent_elev = mean(Tree.pts(cell_ptc_ids,3));
            cell_surf_feaures = super_structure.Component(i).cell(j).surface_features;
            mask = (ismember(cells_plane.cell_ids(:,1), comp_cell_id(1), 'rows'))&...
                   (cell_peaks_features(:,3) >= cell_surf_feaures(3) - threshold.max_distance);
        end
        % Update cell peaks -> 0: inactive
        cells_plane.cell_ids(mask, end) = 0;
     end
end
 
% Re-update inactive cell_peak
mask = ismember(cells_plane.cell_ids(:,[1,2]), inactive_cell_peak_ids(:, [1,2]), 'rows');
cells_plane.cell_ids(mask,3) = 0;
if time_record
    fprintf('Running time of update cells status: %.2f seconds \n', toc);
end