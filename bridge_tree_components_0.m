function [bridge_tree, master_comp_id, source_data] = bridge_tree_components(bridge_tree, tree_level, add_comp_name, master_comp_id, source_comp_id, source_data)
%% This function is to add data points of the components to the data structure of the bridge
% Input:
%     bridge_tree         - Bridge data structure:                             
%     tree_level          - new level of the bridge data structure
%     add_comp_name       - Name of a new component to be added
%     master_comp_id      - Component id that the new one will connect 
%     source_comp_id      - Component ids to be extracted from source and
%     add to the bridge_tree
%     source_data         - Region
% Output:
%     bridge_tree         - Bridge data structure
% Demo: 
% bridge_tree;
% tree_level = 1;
% parent_comp_id = 0;
% child_comp_id = ref_seg_ids;
% source_data = region;
% comp_name = 'Road Surface';


%% Check input parameters
% check level of the bridge tree
if isempty(tree_level) || (tree_level == 0)
    tree_level = length(bridge_tree);    
end

% Preallocation a new level of data structure
if tree_level > length(bridge_tree)
%     current_bridge_tree_level = length(bridge_tree);
    bridge_tree(tree_level).Component = [];
    bridge_tree(tree_level).Link_List = [];
end
% Check parent component ids
if isempty(master_comp_id)
    master_comp_id = 0;
end

% Establish link_list
if isempty(source_comp_id)
    error('Error: Input level must not be empty.')
elseif max(source_comp_id) > length(source_data)
    error('Error: Input child components %d must not higher than the length of source %d.', source_comp_id, length(source_data))
end

if isempty(bridge_tree(tree_level).Link_List)
    new_child_comp_id = 1;
else
%     mask = ismember(bridge_tree(tree_level).Link_List(:,end),slave_comp_id);
%     if any(mask) 
%         new_child_comp_id = unique(bridge_tree(tree_level).Link_List(mask,2));
%     else
        new_child_comp_id = max(bridge_tree(tree_level).Link_List(:,2)) + 1;
%     end
    % Check if parent id is less than new_child_comp_id
%         if parent_comp_id <= new_child_comp_id
%             parent_comp_id = new_child_comp_id;
%             new_child_comp_id = new_child_comp_id + 1;
%         end
end
    new_link_list = zeros(numel(source_comp_id),3);
    new_link_list(:,1) = master_comp_id;
    new_link_list(:,2) = new_child_comp_id;
    new_link_list(:,3) = source_comp_id;
    bridge_tree(tree_level).Link_List = [bridge_tree(tree_level).Link_List;new_link_list];


% clear m n new_link_list

% Add to the structure:
count_cell = 1;
for i = 1:numel(source_comp_id)
    if source_data(source_comp_id(i)).status
        % Update Component and Description
        bridge_tree(tree_level).Component(new_child_comp_id).ids = new_child_comp_id;
        bridge_tree(tree_level).Description(new_child_comp_id).ids = add_comp_name;
        
        for j = 1:length(source_data(source_comp_id(i)).cell)
    %         bridge_tree(bridge_tree_level).Component(new_child_comp_id).ids = child_comp_id(i);
            bridge_tree(tree_level).Component(new_child_comp_id).cell(count_cell).ids = source_data(source_comp_id(i)).cell(j).id;
            bridge_tree(tree_level).Component(new_child_comp_id).cell(count_cell).ptc_ids = source_data(source_comp_id(i)).cell(j).ptc_ids;
            count_cell = count_cell + 1;
        end
        
        % Set the component is transferred to the bridge tree structure
        source_data(source_comp_id(i)).status = 0;
    end
    
end

clear i j count_cell

% Update a new parent id
parent_id = new_child_comp_id;

