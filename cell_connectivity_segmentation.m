function global_region_info = cell_connectivity_segmentation(Tree, cell_ids)
%% This function is to group all cells connection
% Developed by Dr. Linh Truong-Hong, Dept. GRS, TUDelft
% Input:
% Output
%   global_region_info = [
% Demo:
% TREE = OQTR;
% cell_ids = column_cell_info(:,1);
%% Connectivity segmentation
global_region_info = inf(numel(cell_ids),2);
global_region_info(:, 1) = cell_ids;

checked_region_cell = zeros(numel(cell_ids),1);
checked_region_cell(:,1) = cell_ids;

region_no = 0;
while ~isempty(checked_region_cell)
    % Initial seeding voxel
    seed_cell_id = checked_region_cell(1);
    % Seeding region
    seed_region_ids = seed_cell_id;
    seed_region_ids = seed_region_ids(:);
    % Current region
    cur_region = seed_cell_id;
    cur_region = cur_region(:);
    % Update checked voxel
    checked_region_cell(1) = [];
    % Region growing
    while ~isempty(seed_region_ids) 
        
        cur_seed_cell_id = seed_region_ids(1);
        % Searching neighbour
        neighbour_cell_ids = Query_Neighbour_Cells(Tree,cur_seed_cell_id,checked_region_cell);
        
        if ~isempty(neighbour_cell_ids)
            % Update a current region
            cur_region = union(cur_region,neighbour_cell_ids,'rows'); 
            % Remove new add_seeding out of the checked voxel
            mask = ismember(checked_region_cell,neighbour_cell_ids);
            checked_region_cell(mask,:) = [];
            % Add new seeding cell ids
            seed_region_ids = union(seed_region_ids,neighbour_cell_ids,'rows');
        end
        % Remove current_seeding_voxel_ind
        seed_region_ids = setdiff(seed_region_ids,cur_seed_cell_id);
    end
    % Update the current region
    mask = ismember(global_region_info(:,1), cur_region);
    global_region_info(mask,2) = region_no + 1;
    region_no = region_no + 1;
end

%% Update region ids
mask = isinf(global_region_info(:,end));
global_region_info(mask,:) = [];

if ~isempty(global_region_info)
    num_region = max(global_region_info(:,end));
    region_count = histcounts(global_region_info(:,end), num_region);
    region_stat(:, 1) = unique(global_region_info(:,end));
    region_stat(:, 2) = region_count;

    [~, ids] = sort(region_stat(:, 2), "descend");
    region_stat = region_stat(ids,:);

    % Update region_info
    temp_region_ids = inf(size(global_region_info, 1),1);
    for i = 1:size(region_stat,1)
        mask = ismember(global_region_info(:, end), region_stat(i,1));
        temp_region_ids(mask) = i;
    end
    %
    global_region_info(:,end) = temp_region_ids;

    clear num_region region_count region_stat ids mask temp_region_ids
end
