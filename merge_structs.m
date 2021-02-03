function new_data_struct = merge_structs(data_struct,level_ids)
%% This function is to merge data from different level into a new data structure
% Input:
%     level_ids
%     data_struct
% Output:
%     new_data_struct
% Demo:
%     level_ids = components_region_ids
%     data_struct = Region
%% Extract the fields
% There is no function to detect cell id sharing by two segments/regions
new_data_struct = struct('cell',[]);
count = 0;
for i = 1:numel(level_ids)
    level_id = level_ids(i);
    num_cells = length(data_struct(level_id).cell);
    for j = 1:num_cells
        new_data_struct.cell(count + 1).id = data_struct(level_id).cell(j).id;
        new_data_struct.cell(count + 1).ptc_ids = data_struct(level_id).cell(j).ptc_ids;
        new_data_struct.cell(count + 1).surface_features = data_struct(level_id).cell(j).surface_features;
             
        count = count + 1;
    end
    
end

    