function [super_structure, source_data] = bridge_tree_components(super_structure, master_comp_name, slave_comp_id,  slave_comp_name, source_data)
%% This function is to add data points of the components to the data structure of the bridge
% Input:
%     bridge_tree           - Bridge data structure:                             
%     master_comp_id        - Component id that the new one will connect 
%     slave_comp_id         - Component ids to be extracted from source and
%     master_comp_name      - Name of the master component
%     slave_comp_name       - Name of the slave component add to the bridge_tree
%     source_data           - Region
% Output:
%     bridge_tree         - Bridge data structure
%% Demo: 
% supper_structure = SuperStructure;
% master_comp_id = 0;
% slave_comp_id = ref_seg_ids;
% source_data = cells_region;
% add_comp_name = 'Road Surface';


%% Check input parameters
% % Check parent component ids
% if isempty(master_comp_id) ||(master_comp_id ==0)
%     master_comp_name = '';
% end

% Check source_comp_id
if isempty(slave_comp_id)
    error('Error: Input level must not be empty.')
elseif max(slave_comp_id) > length(source_data)
    error('Error: Input child components %d must not higher than the length of source %d.', slave_comp_id, length(source_data))
end

if isempty(super_structure)
    new_comp_id = 1;
else
    new_comp_id = length(super_structure(1).Component) + 1;
end

% Add to the structure:
count_cell = 1;
for i = 1:numel(slave_comp_id)
    if source_data(slave_comp_id(i)).status
        % Update Component and Description
        super_structure(1).Component(new_comp_id).id = new_comp_id;
        super_structure(1).Linked_Component{new_comp_id} = master_comp_name;
        super_structure(1).Description{new_comp_id} = slave_comp_name;

        for j = 1:length(source_data(slave_comp_id(i)).cell)
    %         bridge_tree(bridge_tree_level).Component(new_child_comp_id).ids = child_comp_id(i);
            super_structure(1).Component(new_comp_id).cell(count_cell).ids = source_data(slave_comp_id(i)).cell(j).id;
            super_structure(1).Component(new_comp_id).cell(count_cell).ptc_ids = source_data(slave_comp_id(i)).cell(j).ptc_ids;
            super_structure(1).Component(new_comp_id).cell(count_cell).surface_features = source_data(slave_comp_id(i)).cell(j).surface_features;
            count_cell = count_cell + 1;
        end
        
        % Set the component is transferred to the bridge tree structure
        source_data(slave_comp_id(i)).status = 0;
    end
    
end

clear i j count_cell


