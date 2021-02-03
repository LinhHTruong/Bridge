classdef OctQuadtree < handle
    properties
        pts
        cell_pts
        cell_counts
        cell_bounds
        cell_depths
        cell_parent = zeros(0,1);
        cell_props = ones(1,1);
        tree_prop
    end
    methods
        function this = OctQuadtree(pts,varargin)
%             validateattributes(pts,{'numeric'},{'real','finite','nonnan','ncols', 3})
            validateattributes(pts,{'numeric'},{'real','finite','nonnan'})
            % Allow custom setting of tree_prop
            IP = inputParser;
            IP.addParameter('cell_capacity',0);
            IP.addParameter('max_depth',inf);
            IP.addParameter('max_size',0);
            IP.parse(varargin{:})
            %
            this.tree_prop = IP.Results;
            
            % Initialise a single bin surrounding all given points
            numPts = size(pts,1);
            this.pts = pts;
            this.cell_pts = struct('id',[]);
            this.cell_pts(1,1).id = [1:numPts]';
            this.cell_counts = 1;
            min_bound = min(pts(:,1:3),[],1);
            max_bound = max(pts(:,1:3),[],1);
            bound_length = max_bound - min_bound;
%             this.cell_bounds = [min(pts(:,1:3),[],1) max(pts(:,1:3),[],1)];
%             this.cell_bounds = [min_bound, min_bound+ones(1,3)*max(bound_length)];

                
            if this.tree_prop.max_size ~= 0
                depth = ceil(log2(bound_length(1:2)./this.tree_prop.max_size));
                real_bound_length = (2^max(depth))*this.tree_prop.max_size;
                this.cell_bounds = [min_bound, min_bound+ones(1,3)*real_bound_length];
                this.tree_prop.maxDepth = max(depth);
                max(depth)
                
            else
                this.cell_bounds = [min_bound, min_bound+ones(1,3)*max(bound_length)];
            end
            
            this.cell_depths(1,1) = 0;
            this.cell_parent(1,1) = 0;
            this.cell_props(1,1) = 1;
%             
            %
            this.preallocateSpace;
            this.divide(1);
            this.shrink;
        end
        function preallocateSpace(this)
            numPts = size(this.pts,1);
            if numPts>3
                numBins = 5;
            else
                numBins = 1;
            end
            this.cell_depths(numBins,1) = 0;
            this.cell_parent(numBins,1) = 0;
            this.cell_bounds(numBins,1) = 0; 
        end
        % Recursive Divide
        function divide(this, cell_id)
            for i=1:numel(cell_id)
                cell_no = cell_id(i);
                old_count = this.cell_counts;
                if ~isinf(this.tree_prop.max_depth)
                    % Divide according to a maximum depth criterion
%                     cell_no = cell_id(i);
%                     old_count = this.cell_counts;
                    if (this.cell_depths(cell_no)+1 > this.tree_prop.max_depth)||...
                       (this.cell_props(cell_no)==0)
                        continue;
                    else % Sub-division
                        this.dividecell_(cell_no);
                        this.divide(old_count+1:this.cell_counts);
                        continue;
                    end
                elseif this.tree_prop.max_size ~=0
                    % Divide according to a maximum size of a cell_
                    this_bounds = this.cell_bounds(cell_no,:);
                    cell_edge_size = diff(this_bounds([1:2;4:5])); %Use x and y coordinate
                    minE_edge_size = min(cell_edge_size);
%                     maxEdgeSize = max(cell_edge_size);
                    if (minE_edge_size <= this.tree_prop.max_size)||...
                       (this.cell_props(cell_no)==0)
                        % No subdivision
                        continue;
                    else
                        % Subdivision
                        this.dividecell_(cell_no);
                        this.divide(old_count+1:this.cell_counts);
                        continue;
                    end
                else
                    % Divide according to a maximum number of points within
                    % a cell_
                    if (nnz(this.cell_pts == cell_no) < this.tree_prop.cell_capacity)||...
                       (this.cell_props(cell_no) == 0)
                        % No subdivision
                        continue;
                    else
                        % Subdivision
                        this.dividecell_(cell_no);
                        this.divide(old_count+1:this.cell_counts);
                        continue;
                    end
                end
            end
        end
        % Divide
        function dividecell_(this, cell_no)
            pt_ids = this.cell_pts(cell_no,1).id;            
            cell_pts_xyz = this.pts(pt_ids,:);
            % Get the old corner points and the new division point
            old_min = this.cell_bounds(cell_no, 1:2); %Use x and y coordinates
            old_max = this.cell_bounds(cell_no, 4:5); %Use x and y coordinates
            new_div = mean([old_min; old_max], 1);
            % Build the new boundaries of our 4 subdivisions
            min_mid_max = [old_min, new_div, old_max];
            new_bounds = min_mid_max([...
                                    1 2 3 4;
                                    3 2 5 4;
                                    1 4 3 6;
                                    3 4 5 6]);
            % Determine to which of these 4 bins each current point belongs
            cell_map = cat(3,[0 0],[1 0],...
                             [0 1],[1 1]);
            gt_mask = bsxfun(@gt, cell_pts_xyz(:,1:2), new_div);
            [~,cell_assignment] = max(all(bsxfun(@eq,gt_mask, cell_map), 2),[], 3);
            new_cell_ids = this.cell_counts+1:this.cell_counts + 4;
            % Extract point indicies in new cell_s
            linked_pts_cell = arrayfun(@(i)find(new_cell_ids(cell_assignment)==i), new_cell_ids, 'Un', 0)';
            for j=1:4
                this.cell_pts(new_cell_ids(j),1).id = pt_ids(linked_pts_cell{j});
            end
            this.cell_bounds(new_cell_ids,:) = [new_bounds(:,1:2), repmat(this.cell_bounds(cell_no, 3), [4,1]),...
                                                  new_bounds(:,3:4), repmat(this.cell_bounds(cell_no, 6), [4,1])];
            this.cell_depths(new_cell_ids,1) = this.cell_depths(cell_no) + 1;
            this.cell_parent(new_cell_ids,1) = cell_no;
            % Assigning the property: full = 1; empty = 0;
            new_empty_cell_ids = cellfun(@isempty, linked_pts_cell) == 1;
            this.cell_props(new_cell_ids(new_empty_cell_ids),1) = 0;
            new_full_cell_ids = cellfun(@isempty,linked_pts_cell) == 0;
            this.cell_props(new_cell_ids(new_full_cell_ids),1) = 1;
            this.cell_counts = this.cell_counts + 4;
        end
        %Shrink cell_
        function shrink(this)
            % Shrink height or elevation of the 2D cell_ by z-coordinate of the points within the cell_
            % Extract leaf nodes
            cell_Children = arrayfun(@(i)find(this.cell_parent == i), 1:this.cell_counts, 'Un', 0)';
            leaf_cell_ids = cellfun(@isempty, cell_Children);
%             leaf_cell_ids = leaf_cell_ids(this.cell_props(leaf_cell_ids)==1);
            for i = find(leaf_cell_ids(:))'
                Cell_Shrink_Recurse(i, true)
            end 
            function Cell_Shrink_Recurse(cell_no, mask)
                % Extract the points within the 3D cell_
                old_boundary_min = this.cell_bounds(cell_no,1:3);
                old_boundary_max = this.cell_bounds(cell_no,4:6);
                if mask
                    % Shrink bin based on child POINTS
                    pt_ids = this.cell_pts(cell_no, 1).id;            
                    if this.cell_props(cell_no)==0
                        % Empty cell_, shrink the bin to infinitely small
                        new_bounds = [old_boundary_min, old_boundary_min];
                    else
                        pts_xyz = this.pts(pt_ids,:);
                        new_bounds = [old_boundary_min(:,1:2) max(old_boundary_min(3),min(pts_xyz(:,3))),...
                                      old_boundary_max(:,1:2) min(old_boundary_max(3),max(pts_xyz(:,3)))];
                    end
                else
                    % Shrink bin based on child BINS
                    child_bounds = this.cell_bounds(cell_Children{cell_no},:);
                    new_bounds = [min(child_bounds(:,1:3),[],1), max(child_bounds(:,4:6),[],1)];
                end
                if ~isequal(new_bounds, [old_boundary_min, old_boundary_max])
                    % We just shrunk the boundary. Make it official and
                    % check the parent
                    this.cell_bounds(cell_no,:) = new_bounds;
                    parent_cell_id = this.cell_parent(cell_no);
                    if parent_cell_id>0
                        Cell_Shrink_Recurse(parent_cell_id, false)
                    end
                end
            end
        end
    end
    % cell_ retrieval
    methods
        % Extract leaf nodes
        function leaf_cell_ids = Node_Leaf(this)
            cell_child = arrayfun(@(i)find(this.cell_parent == i), 1:this.cell_counts,'Un',0)';
            leaf_cell_ids = find(cellfun(@isempty, cell_child)==1);
            leaf_cell_ids = leaf_cell_ids(this.cell_props(leaf_cell_ids)==1);
        end
        %Query the point ids within the cell_ Id
        function pts_xyz = queryPoint(this, cell_no)
            if this.cell_props(cell_no) == 0
                pts_xyz = [];
                message('The query cell_ is empty')
            else
                mask = this.cell_pts == cell_no;
                pts_xyz = this.pts(mask,:);
            end
        end
        %Query the cell_id containing input points
        function Cell_Id = Query_Cell_Id(this, cell_ids, query_ptc)
            if size(query_ptc, 1) > 1
                warning('The algorithm is to find the cell_ containing one input point');
                return;
            end
            if isempty(cell_ids)
               cell_ids = Node_Leaf(this); %Retrieve a leaf node cell
            end
            if size(query_ptc, 2) >= 2
                check_cell_bounds = this.cell_bounds(cell_ids,:);
                mask = all(bsxfun(@ge, query_ptc(:,1:2), check_cell_bounds(:,1:2))&...
                           bsxfun(@le, query_ptc(:,1:2), check_cell_bounds(:,4:5)),2);
                cell_ids = cell_ids(mask);
                Cell_Id = unique(cell_ids,'rows');
%                 cell_bounds = this.cell_bounds(query_cell_id,:);
%                 query_ptc = permute(query_ptc(:,1:2),[3 2 1])
%                 [~, cell_NoId] = max(all(bsxfun(@ge, query_ptc,cell_bounds(:,1:2))&...
%                                 bsxfun(@le, query_ptc,cell_bounds(:,4:5)),2),[],1);
%                 cell_NoId = permute(cell_NoId,[1,3,2])';
%                 query_cell_id = query_cell_id(cell_NoId);
%                 cell_Id = unique(query_cell_id,'rows');
            else
                return;
            end
        end
%         % Query a cell_ containing the point
%         function cell_point = querycell_ContainPoint(this,target_cell_ids, pts)
%             if size(pts,1)>1
%                 warning('This function is used to query a cell_ contain one point');
%                 return;
%             else
%                 check_cell_bounds = this.cell_bounds(target_cell_ids,:);
%                 mask = all(bsxfun(@le, pts(:,1:2),check_cell_bounds(:,4:5))&...
%                            bsxfun(@ge, pts(:,4:5),check_cell_bounds(:,1:2)),2);
%                 cell_point = target_cell_ids(mask);
%             end
%         end
        
        function neighbour_cell_ids = Query_Neighbour_Cells(this, query_cell_ids, target_cell_ids)
            flag = false(numel(target_cell_ids),1);
            check_cell_bounds = this.cell_bounds(target_cell_ids,:);
            for i = 1:numel(query_cell_ids)
                query_cell_bounds = this.cell_bounds(query_cell_ids(i),:);
                mask = all(bsxfun(@le, query_cell_bounds(:,1:2), check_cell_bounds(:,4:5))&...
                           bsxfun(@ge, query_cell_bounds(:,4:5), check_cell_bounds(:,1:2)), 2);
                neighbour_cell_ids = target_cell_ids(mask);
                neighbour_cell_ids = neighbour_cell_ids(~ismember(neighbour_cell_ids, query_cell_ids(i)));
                
                mask = ismember(target_cell_ids,neighbour_cell_ids);
                flag(mask) = true;
            end
            neighbour_cell_ids = target_cell_ids(flag);
        end
        
        function neighbour_cell_ids = Query_Neighbour_27Cells(this, query_cell_ids, target_cell_ids)
            flag = false(numel(target_cell_ids),1);
            check_cell_bounds = this.cell_bounds(target_cell_ids,:);
            for i = 1:numel(query_cell_ids)
                query_cell_bounds = this.cell_bounds(query_cell_ids(i),:);
                mask = all(bsxfun(@le, query_cell_bounds(:,1:2), check_cell_bounds(:,4:5))&...
                           bsxfun(@ge, query_cell_bounds(:,4:5), check_cell_bounds(:,1:2)), 2);
                neighbour_cell_ids = target_cell_ids(mask);
                
                mask = ismember(target_cell_ids,neighbour_cell_ids);
                flag(mask) = true;
            end
            neighbour_cell_ids = target_cell_ids(flag);
        end
        
        function neighbour_cell_ids = Query_Edge_Neighbour_Cells(this, query_cell_ids, target_cell_ids)
            flag = false(numel(target_cell_ids),1);
            check_cell_bounds = this.cell_bounds(target_cell_ids,:);
            for i = 1:numel(query_cell_ids)
                query_cell_bounds = this.cell_bounds(query_cell_ids(i),:);
                mask_01 = (bsxfun(@eq, query_cell_bounds(:,1), check_cell_bounds(:,1)))&...
                          (bsxfun(@eq, query_cell_bounds(:,2), check_cell_bounds(:,5))|...
                           bsxfun(@eq, query_cell_bounds(:,5), check_cell_bounds(:,2))); % Vertical
                mask_02 = (bsxfun(@eq, query_cell_bounds(:,2), check_cell_bounds(:,2)))&...
                          (bsxfun(@eq, query_cell_bounds(:,1), check_cell_bounds(:,4))|...
                          bsxfun(@eq, query_cell_bounds(:,4), check_cell_bounds(:,1))); % Horizontal
                mask = mask_01 | mask_02;

                neighbour_cell_ids = target_cell_ids(mask);
                neighbour_cell_ids = neighbour_cell_ids(~ismember(neighbour_cell_ids, query_cell_ids(i)));
                
                mask = ismember(target_cell_ids,neighbour_cell_ids);
                flag(mask) = true;
            end
            neighbour_cell_ids = target_cell_ids(flag);
        end
        
        % Query adjacent cell_s within the window
        function neighbour_cell_ids = Window_Query_Neighbour_Cell(this, query_cell_id, target_cell_ids, window_size)
            if numel(query_cell_id)>1
                warning('This function is used to query a adjacent neighbour of one cell_');
                return;
            else
                check_cell_bounds = this.cell_bounds(target_cell_ids,:);
                query_cell_center = mean([this.cell_bounds(query_cell_id,1:2);
                                        this.cell_bounds(query_cell_id,4:5)],1);
                mask = all(bsxfun(@le, query_cell_center - window_size/2, check_cell_bounds(:,4:5))&...
                           bsxfun(@ge, query_cell_center + window_size/2, check_cell_bounds(:,1:2)), 2);
                neighbour_cell_ids = target_cell_ids(mask);
            end
        end
        % Query adjacent cell_s within a range search
        function neighbour_cell_ids = rangeQuerycell_Neighbour(this, query_cell_id, target_cell_ids,rangesearch)
            if numel(query_cell_id)>1
                warning('This function is used to query a adjacent neighbour of one cell_');
                return;
            else
                check_cell_bounds = this.cell_bounds(target_cell_ids,:);                  
                query_cell_center = mean([this.cell_bounds(query_cell_id,1:2); this.cell_bounds(query_cell_id,4:5)],1);
                mask = all(bsxfun(@le, query_cell_center - rangesearch/2, check_cell_bounds(:,4:5))&...
                           bsxfun(@ge, query_cell_center + rangesearch/2, check_cell_bounds(:,1:2)),2);
                target_cell_ids = target_cell_ids(mask);
                check_cell_center = mean([this.cell_bounds(target_cell_ids,1:2); this.cell_bounds(target_cell_ids,4:5)],1);
                dist_cell_cell_center = sqrt(sum(bsxfun(@power,bsxfun(@minus, query_cell_center,check_cell_center),2),2));                    
                mask = dist_cell_cell_center <= rangesearch;                    
                neighbour_cell_ids = target_cell_ids(mask);
            end
        end
        % Query adjacent lower or equal cell_s of the given cell_
        function neighbour_cell_ids = Query_LE_Neighbour_Cell(this, query_cell_id, target_cell_ids)
            if numel(query_cell_id)>1
                warning('This function is used to query a adjacent neighbour of one cell_');
                return;
            else
                check_cell_bounds = this.cell_bounds(target_cell_ids,:);
                query_cell_bounds = this.cell_bounds(query_cell_id,:);
                mask = all(bsxfun(@le, query_cell_bounds(:,1:2), check_cell_bounds(:,4:5))&...
                           bsxfun(@ge, query_cell_bounds(:,4:5), check_cell_bounds(:,1:2)),2)&...
                       all(bsxfun(@ge, query_cell_bounds(:,6), check_cell_bounds(:,6)),2);
                neighbour_cell_ids = target_cell_ids(mask);
                neighbour_cell_ids = neighbour_cell_ids(~ismember(neighbour_cell_ids, query_cell_id));
            end
        end
        % Query adjacent lower cell_s of the given cell_
        function neighbour_cell_ids = Query_L_Neighbour_Cell(this, query_cell_id, target_cell_ids)
            if numel(query_cell_id)>1
                warning('This function is used to query a adjacent neighbour of one cell_');
                return;
            else
                leaf_cell_ids = Node_Leaf(this); %Retrieve a node cell_
                if ~intersect(query_cell_id,leaf_cell_ids)||~all(intersect(target_cell_ids,leaf_cell_ids))
                    warning('The query cell_ ID is not leaf node. Only leaf node is supported')
                else
                    check_cell_bounds = this.cell_bounds(target_cell_ids,:);
                    query_cell_bounds = this.cell_bounds(query_cell_id,:);
                    mask = all(bsxfun(@le, query_cell_bounds(:,1:2),check_cell_bounds(:,4:5))&...
                               bsxfun(@ge, query_cell_bounds(:,4:5),check_cell_bounds(:,1:2)),2)&...
                           all(bsxfun(@ge, query_cell_bounds(:,3),check_cell_bounds(:,6)),2);
                    neighbour_cell_ids = target_cell_ids(mask);
                    neighbour_cell_ids = neighbour_cell_ids(~ismember(neighbour_cell_ids,query_cell_id));
                end
            end
        end
        % Query adjacent equal cell_s of the given cell_
        function neighbour_cell_ids = Query_Eq_Neighbour_Cell(this, query_cell_id, target_cell_ids)
            if numel(query_cell_id)>1
                warning('This function is used to query a adjacent neighbour of one cell_');
                return;
            else
                leaf_cell_ids = Node_Leaf(this); %Retrieve a node cell_
                if ~intersect(query_cell_id,leaf_cell_ids)||~all(intersect(target_cell_ids, leaf_cell_ids))
                    warning('The query cell_ ID is not leaf node. Only leaf node is supported')
                else
                    check_cell_bounds = this.cell_bounds(target_cell_ids,:);
                    query_cell_bounds = this.cell_bounds(query_cell_id,:);
                    mask = all(bsxfun(@le, query_cell_bounds(:,1:2),check_cell_bounds(:,4:5))&...
                               bsxfun(@ge, query_cell_bounds(:,4:5),check_cell_bounds(:,1:2)),2)&...
                           (all(bsxfun(@ge, query_cell_bounds(:,6),check_cell_bounds(:,6))&...
                                bsxfun(@le, query_cell_bounds(:,3),check_cell_bounds(:,6)),2));
                    neighbour_cell_ids = target_cell_ids(mask);
                    neighbour_cell_ids = neighbour_cell_ids(~ismember(neighbour_cell_ids,query_cell_id));
                end
            end
        end
        % Query equal cell_s of the given cell_. The cell_ having an elevation is the same to one of the query cell_ 
%         function neighbour_cell_ids = QueryEqualcell_(this,query_cell_id,target_cell_ids)
%             if numel(query_cell_id)>1
%                 warning('This function is used to query a adjacent neighbour of one cell_');
%                 return;
%             else
%                 check_cell_bounds = this.cell_bounds(target_cell_ids,:);
%                 query_cell_bounds = this.cell_bounds(query_cell_id,:);
%                 mask = all(bsxfun(@ge, query_cell_bounds(:,6),check_cell_bounds(:,6))&...
%                            bsxfun(@le, query_cell_bounds(:,3),check_cell_bounds(:,6)),2);
%                 neighbour_cell_ids = target_cell_ids(mask);
%                 neighbour_cell_ids = neighbour_cell_ids(~ismember(neighbour_cell_ids,query_cell_id));
%             end
%         end
        %
        % Query the closet cell_ of the given cell_
        %
        function neighbour_cell_ids = Closest_Cell(this, query_cell_id, target_cell_ids)
            if numel(query_cell_id)>1
                warning('This function is used to query a adjacent neighbour of one cell_');
                return;
            else
                check_cell_bounds = this.cell_bounds(target_cell_ids,:);
                check_cell_center = [bsxfun(@plus, check_cell_bounds(:,1), check_cell_bounds(:,4))/2,...
                                 bsxfun(@plus, check_cell_bounds(:,2), check_cell_bounds(:,5))/2,...
                                 bsxfun(@plus, check_cell_bounds(:,3), check_cell_bounds(:,6))/2];
                query_cell_bounds = this.cell_bounds(query_cell_id,:);
                querycell_Center = mean([query_cell_bounds(1:3);query_cell_bounds(4:6)],1);
                d = sqrt(sum(bsxfun(@minus,check_cell_center,querycell_Center).^2,2));
                [~,min_id] = min(d);
                neighbour_cell_ids = target_cell_ids(min_id);
            end
        end
    end
    % Computate feature
    methods
        % The normal OR the largest tangent vector OR the smoothest OR residual  of the cell_
        function cell_features = cell_PCA(this, varargin)
            %feature = 1: normal vector
            %feature = 2: a largest tangent vector
            %feature = 3: The smoothest
            %feature = 4: Residual
            %feature = 5: normal + Residual
            
            if isempty(varargin)
                warning('No plotting since not enough input variables');
                return;
            end
            IP = inputParser;
            IP.addParameter('cell_ids',[]);
            IP.addParameter('num_min_pts',3);
            IP.addParameter('feature',1);
            IP.parse(varargin{:})
            cal_cell_features = IP.Results;
            % Check input
            validateattributes(cal_cell_features.num_min_pts,{'numeric'},{'real','finite','nonnegative','nonnan'})
            if (cal_cell_features.feature<1)||(cal_cell_features.feature>5)
                cal_cell_features.feature = 1; %Cal normal vector
            end
            % Calculate feature
            ref_vector = [0,0,1];
            if isempty(cal_cell_features.cell_ids)
                cell_ids = Node_Leaf(this);
            else
                cell_ids = cal_cell_features.cell_ids;
            end
            
            % Prelocation 
            if (cal_cell_features.feature == 1) ||(cal_cell_features.feature == 2)
                cell_vector = zeros(numel(cell_ids),4);
            elseif  (cal_cell_features.feature == 3) ||(cal_cell_features.feature == 4)
                cell_vector = zeros(numel(cell_ids),2);
            else
                cell_vector = zeros(numel(cell_ids),5);
            end
            cell_vector(:,1) = cell_ids;
            for i=1:numel(cell_ids)
                mask = this.cell_pts(cell_ids(i)).id;
                cell_pts_xyz = this.pts(mask,:);
                if size(cell_pts_xyz,1) < cal_cell_features.num_min_pts
                    cell_vector(i,2:4) = inf;
                else
                    if cal_cell_features.feature == 1 %Normal vector
                        normalvector = eigenspace(cell_pts_xyz,1);
                        if sign(dot(normalvector,ref_vector)) > 0
                            cell_vector(i,2:4) = normalvector;
                        else
                            cell_vector(i,2:4) = -normalvector;
                        end
                    elseif cal_cell_features.feature == 2 %The tangent vector
                        cell_vector(i,2:4) = eigenspace(cell_pts_xyz,2);
                    elseif cal_cell_features.feature == 3 %The smoothest
                        cell_vector(i,2) = eigenspace(cell_pts_xyz,3);
                    elseif  cal_cell_features.feature == 4 %Residual
                        normVect = eigenspace(cell_pts_xyz,1);
                        res = cell_pts_xyz*(normVect)'+(-mean(cell_pts_xyz)*(normVect)');
                        cell_vector(i,2)= sqrt((1/size(res,1))*sum(res.^2));
                    else % cal_cell_features.feature == 4 : normal + Residual
                        normVect = eigenspace(cell_pts_xyz,1);
                        res = cell_pts_xyz*(normVect)'+(-mean(cell_pts_xyz)*(normVect)');
                        cell_vector(i,2:4) = normVect;
                        cell_vector(i,5)= sqrt((1/size(res,1))*sum(res.^2));
                    end
                end
            end
            mask = isinf(cell_vector(:,2));
            cell_features = cell_vector(~mask,:);
        end
        % The Height diff Or Z gradient Or both
        function cell_elev = Call_Cell_zFeatures(this, varargin)
            % Z different = 1
            % Z Gradient = 2
            % Slope = 3
            % Both = 4
            IP = inputParser;
            IP.addParameter('feature',1);
            IP.parse(varargin{:})
            cell_zfeatures = IP.Results;
            if (cell_zfeatures.feature < 0)||(cell_zfeatures.feature >4)
                cell_zfeatures.feature = 1;
            end
            % Calculate feature
            leaf_cell_ids = Node_Leaf(this);
            cell_z = zeros(numel(leaf_cell_ids), 4);
            cell_z(:,1) = leaf_cell_ids;
            for i=1:numel(leaf_cell_ids)
                mask = this.cell_pts == leaf_cell_ids(i);
                cell_pts_xyz = this.pts(mask,:);
                if size(cell_pts_xyz,1) <= 1
                    cell_z(i,2:3) = inf;
                else
                    if cell_zfeatures.feature == 1
                        cell_z(i,2) = diff([min(cell_pts_xyz(:,3)),max(cell_pts_xyz(:,3))]);   
                    elseif cell_zfeatures.feature == 2
                        [~,id_max] = max(cell_pts_xyz(:,3));
                        [~,id_min] = min(cell_pts_xyz(:,3));
                        cell_z(i,2) = diff([cell_pts_xyz(id_min,3),cell_pts_xyz(id_max,3)])/...
                                     norm(cell_pts_xyz(id_max,1:2)-cell_pts_xyz(id_min,1:2));
                    else
                        [~,id_max] = max(cell_pts_xyz(:,3));
                        [~,id_min] = min(cell_pts_xyz(:,3));
                        cell_z(i,2) = diff([cell_pts_xyz(id_min,3),cell_pts_xyz(id_max,3)]);
                        cell_z(i,3) = diff([cell_pts_xyz(id_min,3),cell_pts_xyz(id_max,3)])/...
                                     norm(cell_pts_xyz(id_max,1:2)-cell_pts_xyz(id_min,1:2));
                        ref_vector = [0,0,1];
                        normalvector = eigenspace(cell_pts_xyz,1);
                        if sign(dot(normalvector,ref_vector)) < 0
                            normalvector = -normalvector;
                        end
                        cell_z(i,4) = tan(dot(normalvector,ref_vector));          
                    end
                end
            end
            cell_elev = cell_z(~isinf(cell_z(:,2)),:);
            if cell_zfeatures.feature ~= 4
                cell_elev(:,3:4) = [];
            end
        end
    end
    %Visualization
    methods %Visualization
        function plotAll(this, varargin)
            if isempty(varargin)
                warning('No plotting since not enough input variables');
                return;
            end
            IP = inputParser;
            IP.addParameter('cell_Depth',[]);
            IP.addParameter('cell_Id',[]);
            IP.addParameter('typePlotcell_',2);
            IP.addParameter('typePlotPoint',0);
            IP.addParameter('color',[]);
            IP.addParameter('fill',1); % 1: filled cubic; 0; wireframe
            IP.parse(varargin{:})
            plotOption = IP.Results;
            %
            if isempty(plotOption.cell_Depth)&&isempty(plotOption.cell_Id)
                cell_Id = Node_Leaf(this);
            elseif ~isempty(plotOption.cell_Depth)&& isempty(plotOption.cell_Id)
                if (plotOption.cell_Depth<0)||(plotOption.cell_Depth>max(this.cell_depths))
                    plotOption.cell_Depth = max(this.cell_depths);
                end
                cell_Id = this.cell_depths == plotOption.cell_Depth;
            else
                cell_Id = plotOption.cell_Id;
            end
            %
            if (plotOption.typePlotcell_<0)||(plotOption.typePlotcell_>2)
                typePlotcell_ = 2; %Both empty and full cell_
            else
                typePlotcell_ = plotOption.typePlotcell_; %Either empty "0" or full "1"
            end
            %
            if (plotOption.typePlotPoint<0)||(plotOption.typePlotPoint>1)
                typePlotPoint = 0; %Default no plot points
            else
                typePlotPoint = plotOption.typePlotPoint; %Either plot points "1" or not "0"
            end
            %
            if ischar(plotOption.color)
                color = rem(floor((strfind('kbgcrmyw', plotOption.color) - 1)*[0.25 0.5 1]), 2);
            else
                color = plotOption.color;
            end
            %
            if (plotOption.fill<0)||(plotOption.fill>1)
                fill = 1;
            else
                fill = plotOption.fill;
            end
            %
            if typePlotcell_~=2 %Either empty "0" or full "1"
                cell_Id = cell_Id(this.cell_props(cell_Id)==typePlotcell_);
            end
            listcell_ = (1:this.cell_counts)'; 
            cell_Id = listcell_(cell_Id);
            % Plot
            plotcell_(this, cell_Id, typePlotcell_, color, fill)
            if (typePlotcell_~=0)&&(typePlotPoint==1)&& (fill==0)
                plotPoint(this,cell_Id, color) 
            end
        end
        % Plot cell_   
        function plotcell_(this, cell_Id, typePlotcell_, color, fill)
            hold on;
            if isempty(color)
               color = lines(numel(cell_Id));
            else
                color = repmat(color,[numel(cell_Id),1]);
            end
            for i = 1:numel(cell_Id)
                cell_MinMax = this.cell_bounds(cell_Id(i),:);
                if fill==1 % Plot filled cell_   
                    face1 = cat(1,cell_MinMax([1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3]));
                    face2 = cat(1,cell_MinMax([1 2 6; 1 5 6; 4 5 6; 4 2 6; 1 2 6]));
                    face3 = cat(1,cell_MinMax([1 2 3; 1 2 6; 4 2 6; 4 2 3; 1 2 3]));
                    face4 = cat(1,cell_MinMax([4 2 3; 4 2 6; 4 5 6; 4 5 3; 4 2 3]));
                    face5 = cat(1,cell_MinMax([4 5 3; 4 5 6; 1 5 6; 1 5 3; 4 5 3]));
                    face6 = cat(1,cell_MinMax([1 5 3; 1 5 6; 1 2 6; 1 2 3; 1 5 3]));
                    if typePlotcell_==2
                        if this.cell_props(cell_Id(i))==1
                            fill3(face1(:,1),face1(:,2),face1(:,3),color(i,:))
                            fill3(face2(:,1),face2(:,2),face2(:,3),color(i,:))
                            fill3(face3(:,1),face3(:,2),face3(:,3),color(i,:))
                            fill3(face4(:,1),face4(:,2),face4(:,3),color(i,:))
                            fill3(face5(:,1),face5(:,2),face5(:,3),color(i,:))
                            fill3(face6(:,1),face6(:,2),face6(:,3),color(i,:))
                        else
                            fill3(face1(:,1),face1(:,2),face1(:,3),'w')
                            fill3(face2(:,1),face2(:,2),face2(:,3),'w')
                            fill3(face3(:,1),face3(:,2),face3(:,3),'w')
                            fill3(face4(:,1),face4(:,2),face4(:,3),'w')
                            fill3(face5(:,1),face5(:,2),face5(:,3),'w')
                            fill3(face6(:,1),face6(:,2),face6(:,3),'w')
                        end
                    else
                        fill3(face1(:,1),face1(:,2),face1(:,3),color(i,:))
                        fill3(face2(:,1),face2(:,2),face2(:,3),color(i,:))
                        fill3(face3(:,1),face3(:,2),face3(:,3),color(i,:))
                        fill3(face4(:,1),face4(:,2),face4(:,3),color(i,:))
                        fill3(face5(:,1),face5(:,2),face5(:,3),color(i,:))
                        fill3(face6(:,1),face6(:,2),face6(:,3),color(i,:))
                    end   
                else %No fill cell_
                    cell_Edge = cat(1, cell_MinMax([...
                        1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3;...
                        1 2 6; 4 2 6; 4 5 6; 1 5 6; 1 2 6; 1 2 3]),...
                        nan(1,3), cell_MinMax([4 2 3; 4 2 6]),...
                        nan(1,3), cell_MinMax([4 5 3; 4 5 6]),...
                        nan(1,3), cell_MinMax([1 5 3; 1 5 6]));
                    if typePlotcell_==2
                        if this.cell_props(cell_Id(i))==1
                           plot3(cell_Edge(:,1),cell_Edge(:,2),cell_Edge(:,3),...
                               'Color',color(i,:),'LineWidth', 2);
                        else
                            plot3(cell_Edge(:,1),cell_Edge(:,2),cell_Edge(:,3),...
                               'Color','k','LineWidth', 2);
                        end
                    else
                        plot3(cell_Edge(:,1),cell_Edge(:,2),cell_Edge(:,3),...
                               'Color',color(i,:),'LineWidth', 2);
                    end
                end
            end
        end 
        % Plot points
        function plotPoint(this,cell_Id, color)
            hold on
            doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
            if isempty(color)
               color = lines(numel(cell_Id));
            else
                color = repmat(color,[numel(cell_Id),1]);
            end
            %
            for i = 1:numel(cell_Id)
                
                doplot3(this.pts(this.cell_pts(cell_Id(i)).id,:),'.','Color',color(i,:))
            end 
        end
 
        function plotOrdercell_(this)
            max_depth = max(this.cell_depths);
            numcell_ = this.cell_counts;
            cell_Id = 1:numcell_;
            h = zeros(numcell_,1);
            count = 0;
            hold all;
            for i=1:max_depth+1
                mask = this.cell_depths==(i-1);
                cell_Id_level = cell_Id(mask);
                colors = lines(numel(cell_Id_level));
                for j=1:numel(cell_Id_level)
                    count = count+1;
                    cell_MinMax = this.cell_bounds(cell_Id_level(j),:);
                    cell_Edge = cat(1, cell_MinMax([...
                                    1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3;...
                                    1 2 6; 4 2 6; 4 5 6; 1 5 6; 1 2 6; 1 2 3]),...
                                    nan(1,3), cell_MinMax([4 2 3; 4 2 6]),...
                                    nan(1,3), cell_MinMax([4 5 3; 4 5 6]),...
                                    nan(1,3), cell_MinMax([1 5 3; 1 5 6]));
                    h(count) = plot3(cell_Edge(:,1),cell_Edge(:,2),cell_Edge(:,3));
                    set(h(count),'Color',colors(j,:),'LineWidth', 2)
                    pause
                end
                
            end
        end
        function plotcell_Attribute(this, varargin)
            IP = inputParser;
            IP.addParameter('cell_Id',[]);
            IP.addParameter('attribute',[]);
            IP.addParameter('colormap','hsv');
            IP.addParameter('fill',1); % 1: filled cubic; 0; wireframe
            IP.parse(varargin{:})
            plotOption = IP.Results;
            if isempty(plotOption.cell_Id)||isempty(plotOption.attribute)
                return;
            end
            %
            cell_Id = plotOption.cell_Id;
            attribute = plotOption.attribute;
            cmap = colormap(plotOption.colormap);
            normalizeAtt = (attribute-min(attribute))/(max(attribute)-min(attribute));
            rangeAtt = linspace(min(normalizeAtt),max(normalizeAtt),size(cmap,1)+1);
            for i=1:numel(rangeAtt)-1
                mask = (normalizeAtt>=rangeAtt(i))&(normalizeAtt<rangeAtt(i+1));
                curcell_Id = cell_Id(mask);
                if ~isempty(curcell_Id)
                    for j=1:numel(curcell_Id)
                        cell_MinMax = this.cell_bounds(curcell_Id(j),:);
                        if plotOption.fill==1
                            face1 = cat(1,cell_MinMax([1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3]));
                            face2 = cat(1,cell_MinMax([1 2 6; 1 5 6; 4 5 6; 4 2 6; 1 2 6]));
                            face3 = cat(1,cell_MinMax([1 2 3; 1 2 6; 4 2 6; 4 2 3; 1 2 3]));
                            face4 = cat(1,cell_MinMax([4 2 3; 4 2 6; 4 5 6; 4 5 3; 4 2 3]));
                            face5 = cat(1,cell_MinMax([4 5 3; 4 5 6; 1 5 6; 1 5 3; 4 5 3]));
                            face6 = cat(1,cell_MinMax([1 5 3; 1 5 6; 1 2 6; 1 2 3; 1 5 3]));
                            %
                            fill3(face1(:,1),face1(:,2),face1(:,3),cmap(i,:))
                            fill3(face2(:,1),face2(:,2),face2(:,3),cmap(i,:))
                            fill3(face3(:,1),face3(:,2),face3(:,3),cmap(i,:))
                            fill3(face4(:,1),face4(:,2),face4(:,3),cmap(i,:))
                            fill3(face5(:,1),face5(:,2),face5(:,3),cmap(i,:))
                            fill3(face6(:,1),face6(:,2),face6(:,3),cmap(i,:))
                        else
                            cell_Edge = cat(1, cell_MinMax([...
                                1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3;...
                                1 2 6; 4 2 6; 4 5 6; 1 5 6; 1 2 6; 1 2 3]),...
                                nan(1,3), cell_MinMax([4 2 3; 4 2 6]),...
                                nan(1,3), cell_MinMax([4 5 3; 4 5 6]),...
                                nan(1,3), cell_MinMax([1 5 3; 1 5 6]));
                            plot3(cell_Edge(:,1),cell_Edge(:,2),cell_Edge(:,3),...
                               'Color',cmap(i,:),'LineWidth', 2);
                        end
                    end         
                end
            end
%             h = colorbar;
            cbh = colorbar('YGrid','on');
            num_space = size(get(cbh,'ytick'),2)-1;
            colorLabel = cell_(1,num_space+1);
            plot_dist = round(linspace(min(attribute),max(attribute),num_space+1)*100)/100;
            for j=1:numel(plot_dist)
                colorLabel{j} = sprintf('%0.2f',plot_dist(j));
            end
            set(cbh,'yticklabel',colorLabel,'FontSize',10,'FontName','Time New Roman');
        end
    end %End a method visualization
end