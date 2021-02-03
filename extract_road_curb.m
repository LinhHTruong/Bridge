function super_structure = extract_road_curb(Tree, super_structure, threshold)
% This function is to extract the road curb located between the road
% surface and footpath
% Input:
%     Tree                - Data structure stored 2D cells
%     bridge_tree         - Data structure of the bridge stored a final result
% Output:
%     bridge_tree

% Demo:

% Tree
% super_structure  = SuperStructure
% cells_region;

%% Extract indices of the road and foot path
% Road surface
mask = cellfun(@(s) contains('Road Surface', s), super_structure(1).Description);
road_surf_id = find(mask);
% Foot path
mask = cellfun(@(s) contains(s, 'Footpath'), super_structure(1).Description);
footpath_ids = find(mask);
%% Extract the road curb
if ~isempty(road_surf_id)
    if ~isempty(footpath_ids)
        
        % Retrieve road surface
        road_cell_ids = vertcat(super_structure(1).Component(road_surf_id).cell.ids);
        road_cell_surf = vertcat(super_structure(1).Component(road_surf_id).cell.surface_features);
        [road_bound_cell_ids, ~] = get_boundary_cells(Tree, road_cell_ids(:,1));
        mask = ismember(road_cell_ids(:,1), road_bound_cell_ids);
        road_bound_cell_surf = road_cell_surf(mask,:);
        clear road_cell_surf

        % Check for each footpath
        for i = 1: numel(footpath_ids)

            % Retrieve segments connected to the ref_seg
            footpath_id = footpath_ids(i);
            footpath_cell_ids = vertcat(super_structure(1).Component(footpath_id).cell.ids);
            footpath_cell_surf = vertcat(super_structure(1).Component(footpath_id).cell.surface_features);
            [footpath_bound_cell_ids, ~] = get_boundary_cells(Tree, footpath_cell_ids(:,1));        
            mask = ismember(footpath_cell_ids(:,1), footpath_bound_cell_ids);
            footpath_bound_cell_surf = footpath_cell_surf(mask,:);
            clear footpath_cell_surf

            % Searching cells can contain points of the curb
            % For road
            road_curb_cells_ids = Query_Edge_Neighbour_Cells(Tree, footpath_bound_cell_ids, road_bound_cell_ids);
            mask = ismember(road_bound_cell_ids, road_curb_cells_ids);
            road_curb_cells_surf = road_bound_cell_surf(mask,:);
            
            % For Footpath
            footpath_curbs_cell_ids = Query_Edge_Neighbour_Cells(Tree, road_curb_cells_ids, footpath_bound_cell_ids);
            mask = ismember(footpath_bound_cell_ids, footpath_curbs_cell_ids);
            footpath_curbs_cells_surf = footpath_bound_cell_surf(mask,:);
            
            % Assembly
            curb_cells_ids = union(road_curb_cells_ids, footpath_curbs_cell_ids, 'stable');
            curb_cells_surf = union(road_curb_cells_surf, footpath_curbs_cells_surf, 'rows', 'stable');

            if isempty(curb_cells_ids)
                continue;
            end
            %% Check for each cell to find a road curb in one side
            % Searching the road curb
            road_curb = struct('cell',[],'status',[]);
            road_curb.status = 1;
            no_curb_cell = 1;
            for j = 1:size(curb_cells_ids,1)
                
                % Retrieve cells
                curb_cell_id = curb_cells_ids(j);
                
                % Classify cells: 1)Contain points of road and footpath; 2)
                % only road; 3 only footpath; Road is considered as the
                % first surface while the footpath is the second
                if find(ismember(road_curb_cells_ids, curb_cell_id) ==1)
                    road_curb_cells_id = curb_cell_id;
                else
                    neighbour_cell_local_id = knnsearch(road_curb_cells_surf(:,1:3), curb_cells_surf(i,1:3), 'k', 1);
                    road_curb_cells_id = road_curb_cells_ids(neighbour_cell_local_id);
                    clear neighbour_cell_local_id
                end
                % Determine the first surface based the road surface
                mask = ismember(road_cell_ids(:,1), road_curb_cells_id);
                st_cell_surf = super_structure(1).Component(road_surf_id).cell(mask).surface_features;
                st_cell_pts_ids = super_structure(1).Component(road_surf_id).cell(mask).ptc_ids;

                
                % For the footpath
                if find(ismember(footpath_curbs_cell_ids, curb_cell_id) ==1)
                    footpath_curb_cells_id = curb_cell_id;
                else
                    neighbour_cell_local_id = knnsearch(footpath_curbs_cells_surf(:,1:3), curb_cells_surf(i,1:3), 'k', 1);
                    footpath_curb_cells_id = footpath_curbs_cell_ids(neighbour_cell_local_id);
                    clear neighbour_cell_local_id
                    
                end
                % Determine the first surface based the footpath surface
                mask = ismember(footpath_cell_ids(:,1), footpath_curb_cells_id);
                nd_cell_surf = super_structure(1).Component(footpath_id).cell(mask).surface_features;
                nd_cell_pts_ids = super_structure(1).Component(footpath_id).cell(mask).ptc_ids;
                
                % Extract the points between two surfaces: road and footpath
                cell_pts_ids = Tree.cell_pts(curb_cell_id).id;
                cell_pts_xyz = Tree.pts(cell_pts_ids,1:3);

                % Extract the points between the surfaces
                [ptc_local_ids, ~] = extract_ptc_2_planes(cell_pts_xyz, st_cell_surf, nd_cell_surf);
                road_curb_cell_pts_ids = cell_pts_ids(ptc_local_ids);
                if numel(road_curb_cell_pts_ids) < threshold.min_num_pts
                    continue;
                end
                
                % Remove points within two surface
                road_curb_cell_pts_ids = setdiff(road_curb_cell_pts_ids, union(st_cell_pts_ids, nd_cell_pts_ids));
                if numel(road_curb_cell_pts_ids) < threshold.min_num_pts
                    continue;
                end    
                % Filter the points on the road curb
                mask = ismember(cell_pts_ids, road_curb_cell_pts_ids);
                road_curb_cell_pts_xyz = cell_pts_xyz(mask,:);
                pts_local_id = road_curb_ransac_filter(road_curb_cell_pts_xyz, threshold);
                if numel(pts_local_id) < threshold.min_num_pts
                    continue;
                end  
                road_curb_cell_pts_ids = road_curb_cell_pts_ids(pts_local_id);
                road_curb_cell_pts_xyz = road_curb_cell_pts_xyz(pts_local_id,1:3);
                % Update road curb
                if numel(road_curb_cell_pts_ids) >= threshold.min_num_pts
                    road_curb.cell(no_curb_cell).id = [curb_cell_id, 0];
                    road_curb.cell(no_curb_cell).ptc_ids = road_curb_cell_pts_ids;
                    road_curb.cell(no_curb_cell).surface_features = [mean(road_curb_cell_pts_xyz), eigenspace(road_curb_cell_pts_xyz,1)]; 
                    no_curb_cell = no_curb_cell + 1;
                end
            end
            
            %% Update the bridge tree
            master_comp_name = super_structure(1).Description{road_surf_id};
            footpath_name = super_structure(1).Description{footpath_id};
            footpath_num = sscanf(footpath_name,'Footpath%d');
            slave_comp_id = 1;
            slave_comp_name = strcat("Road Curb ", num2str(footpath_num));
            [super_structure, ~] = bridge_tree_components(super_structure, master_comp_name, slave_comp_id, slave_comp_name, road_curb);


        end
    end
end
        