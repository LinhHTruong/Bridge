classdef OcTree < handle
    properties
        pts
        voxel_ptc_ids
        voxel_count
        voxel_bounds
        depth
%         voxel_size
        parents = zeros(0,1);
        voxel_prop = ones(1,1);
        octree_property
    end
    methods %Establish Octree
        function this = OcTree(pts,varargin)
            validateattributes(pts,{'numeric'},{'real','finite','nonnan','ncols', 3})
            % Initialise a single bin surrounding all given points
            numPts = size(pts,1);
            this.pts = pts;
%             this.voxel_ptc_ids = ones(numPts,1);
            this.voxel_ptc_ids = struct('id',[]);
            this.voxel_ptc_ids(1,1).id = [1:numPts]';
            this.voxel_count = 1;
            min_bound = min(pts(:,1:3),[],1);
            max_bound = max(pts(:,1:3),[],1);
            bound_length = max_bound - min_bound;
            this.voxel_bounds = [min_bound, min_bound+ones(1,3)*max(bound_length)];
            
%             this.voxel_bounds = [min(pts,[],1) max(pts,[],1)];
            this.depth(1,1) = 0;
%             this.voxel_size = 0;
            this.parents(1,1) = 0;
            this.voxel_prop(1,1) = 1;
%             % Allow custom setting of octree_property
            IP = inputParser;
            IP.addParameter('voxelCapacity',0);
            IP.addParameter('maxDepth',inf);
            IP.addParameter('maxSize',0);
            IP.addParameter('minNoPts',5);
            IP.addParameter('Decompose',1);
            IP.addParameter('min_ptc',5);
            IP.addParameter('normal_direction', 'all');
            IP.parse(varargin{:})
            this.octree_property = IP.Results;
            this.preallocateSpace;
            this.divide(1);
        end
        %
        function preallocateSpace(this)
            numPts = size(this.pts,1);
            if numPts>3
                numBins = 9;
            else
                numBins = 1;
            end
            this.depth(numBins,1) = 0;
            this.parents(numBins,1) = 0;
            this.voxel_bounds(numBins,1) = 0; 
        end
        % Divide
        function divide(this, startingVoxel)
            for i=1:numel(startingVoxel)
                voxelNo = startingVoxel(i);
                oldCount = this.voxel_count;
                if ~isinf(this.octree_property.maxDepth)
                    % Divide according to a maximum voxel_bounds criterion
%                     voxelNo = startingVoxel(i);
%                     oldCount = this.voxel_count;
                    if (this.depth(voxelNo)+1 > this.octree_property.maxDepth)||...
                       (this.voxel_prop(voxelNo)==0)
                        continue;
                    else % Sub-division
                        this.divideVoxel(voxelNo);
                        this.divide(oldCount+1:this.voxel_count);
                        continue;
                    end
                elseif this.octree_property.maxSize ~=0
                    % Divide according to a maximum size of a voxel
                    thisBounds = this.voxel_bounds(voxelNo,:);
                    voxelEdgeSize = diff(thisBounds([1:3;4:6]));
                    if strcmp(this.octree_property.normal_direction, 'xy')||strcmp(this.octree_property.normal_direction, 'yx')
                        cell_edge_size = max(voxelEdgeSize([1,2])); %Use x and y coordinate
                    elseif strcmp(this.octree_property.normal_direction, 'yz')||strcmp(this.octree_property.normal_direction, 'zy')
                        cell_edge_size = max(voxelEdgeSize([2,3])); %Use y and z coordinate
                    elseif strcmp(this.octree_property.normal_direction, 'xz')||strcmp(this.octree_property.normal_direction, 'zx')
                        cell_edge_size = max(voxelEdgeSize([1,3])); %Use y and z coordinate
                    else
                        cell_edge_size = max(voxelEdgeSize); %Use x and z coordinate
                    end

                    maxE_edge_size = max(cell_edge_size);
                    
%                     EdgeSize = min(voxelEdgeSize);
%                     EdgeSize = max(voxelEdgeSize);
                    if (maxE_edge_size < this.octree_property.maxSize)||...
                       (this.voxel_prop(voxelNo)==0)
                        % No subdivision
                        continue;
                    else
                        % Subdivision
                        this.divideVoxel(voxelNo);
                        this.divide(oldCount+1:this.voxel_count);
                        continue;
                    end
                else
                    % Divide according to a maximum number of points within
                    % a voxel
                    if (nnz(this.voxel_ptc_ids==voxelNo) < this.octree_property.voxelCapacity)||...
                       (this.voxel_prop(voxelNo)==0)
                        % No subdivision
                        continue;
                    else
                        % Subdivision
                        this.divideVoxel(voxelNo);
                        this.divide(oldCount+1:this.voxel_count);
                        continue;
                    end
                end
            end
        end
        %
        function divideVoxel(this, voxelNo)
            ptsId = this.voxel_ptc_ids(voxelNo,1).id;            
            voxelPts = this.pts(ptsId,:);
            % Get the old corner points and the new division point
            oldMin = this.voxel_bounds(voxelNo,1:3);
            oldMax = this.voxel_bounds(voxelNo,4:6);
            newDiv = mean([oldMin; oldMax], 1);
            % Build the new boundaries of our 8 subdivisions
            minMidMax = [oldMin newDiv oldMax];
            newBounds = minMidMax([...
                                    1 2 3 4 5 6;
                                    1 2 6 4 5 9;
                                    1 5 3 4 8 6;
                                    1 5 6 4 8 9;
                                    4 2 3 7 5 6;
                                    4 2 6 7 5 9;
                                    4 5 3 7 8 6;
                                    4 5 6 7 8 9]);
            % Determine to which of these 8 bins each current point belongs
            voxelMap = cat(3,[0 0 0],[0 0 1],[0 1 0],[0 1 1],...
                             [1 0 0],[1 0 1],[1 1 0],[1 1 1]);
            gtMask = bsxfun(@gt, voxelPts, newDiv);
            [~,voxelAssignment] = max(all(bsxfun(@eq,gtMask,voxelMap),2),[],3);
            newvoxelId = [this.voxel_count+1:this.voxel_count+8]';
            %Extract point indicies in new voxels
            linkedPtsVoxel = arrayfun(@(i)find(newvoxelId(voxelAssignment)==i),newvoxelId,'Un',0)';
%             [this.voxel_ptc_ids(newvoxelId,1).id] = deal(linkedPtsVoxel{:});
            for j=1:8
                this.voxel_ptc_ids(newvoxelId(j),1).id = ptsId(linkedPtsVoxel{j});
            end
            this.voxel_bounds(newvoxelId,:) = newBounds;
            this.depth(newvoxelId,1) = this.depth(voxelNo)+1;
            this.parents(newvoxelId,1) = voxelNo;
            % Assigning the property: full = 1; empty = 0;
%             newEmptyVoxel = cellfun(@isempty,linkedPtsVoxel)>=0;
            newEmptyVoxel = cellfun(@length,linkedPtsVoxel) < this.octree_property.min_ptc;
            this.voxel_prop(newvoxelId(newEmptyVoxel),1) = 0;
%             newFullVoxel = cellfun(@isempty,linkedPtsVoxel)==0;
            newFullVoxel = cellfun(@length,linkedPtsVoxel) >= this.octree_property.min_ptc;
            this.voxel_prop(newvoxelId(newFullVoxel),1) = 1;
            this.voxel_count = this.voxel_count + 8;
        end
    end
    methods %Retrive data
        % Extract leaf nodes -> full voxels
        function voxelIsLeaf = leafNode(this)
            voxelChildren = arrayfun(@(i)find(this.parents==i),1:this.voxel_count,'Un',0)';
            voxelIsLeaf = find(cellfun(@isempty,voxelChildren)==1);
            voxelIsLeaf = voxelIsLeaf(this.voxel_prop(voxelIsLeaf)==1);
        end
        % Extract leaf nodes -> empty voxels
        function voxelIsLeaf = leafNodeEmpty(this)
            % Leaf node is no childrent
            voxelChildren = arrayfun(@(i)find(this.parents==i),1:this.voxel_count,'Un',0)';
            voxelIsLeaf = find(cellfun(@isempty,voxelChildren)==1);
            % Leaf node: empty voxel is prop == 0
            voxelIsLeaf = voxelIsLeaf(this.voxel_prop(voxelIsLeaf)==0);
            % Leaf node is empty voxels contain point(s) -> lenght > 0
            vocelIsLeaf_ptc = (arrayfun(@(x) vertcat(this.voxel_ptc_ids(x).id),voxelIsLeaf,'UniformOutput',false));
            mask = cellfun(@length,vocelIsLeaf_ptc) > 0;
            voxelIsLeaf = voxelIsLeaf(mask);
        end
        %Query the voxel id containing input points
        function voxelId = queryVoxel(this, newPts)
            if size(newPts,2)==3
                voxelLeaf = leafNode(this); %Retrieve a node voxel
                voxelLeafBoundaries = this.voxel_bounds(voxelLeaf,:);
                newPts = permute(newPts,[3 2 1]);
                [~, voxelNoId] = all(bsxfun(@ge, newPts,voxelLeafBoundaries(:,1:3))&...
                                     bsxfun(@le, newPts,voxelLeafBoundaries(:,4:6)),2);
                voxelNoId = permute(voxelNoId,[1,3,2])';
                voxelId = voxelLeaf(voxelNoId);
            else
                return;
            end
        end
        % Query 26 adjacent voxels of the given voxel
        function neighbour_voxel_ids = query26_neighbour_voxels(this, query_voxel_ids, target_voxel_ids)
            
            flag = false(numel(target_voxel_ids),1);
            target_voxel_bounds = this.voxel_bounds(target_voxel_ids,:);
            for i = 1:numel(query_voxel_ids)
                query_voxel_bounds = this.voxel_bounds(query_voxel_ids(i),:);
                mask = all(bsxfun(@le, query_voxel_bounds(:,1:3),target_voxel_bounds(:,4:6))&...
                           bsxfun(@ge, query_voxel_bounds(:,4:6),target_voxel_bounds(:,1:3)),2);
                neighbour_voxel_ids = target_voxel_ids(mask);
                neighbour_voxel_ids = setdiff(neighbour_voxel_ids, query_voxel_ids(i));
                
                mask = ismember(target_voxel_ids, neighbour_voxel_ids);
                flag(mask) = true;
            end
            neighbour_voxel_ids = target_voxel_ids(flag);
        end
        % Query 26 adjacent voxels of the given voxel
        function voxelNeighbourId = query27VoxelNeighbour(this, queryVoxelId, sourceVoxelId)
            if (numel(queryVoxelId)>1) || isempty(sourceVoxelId)
                warning('This function is used to query a adjacent neighbour of one voxel or target voxel should not empty');
                return;
            else
%                 if isempty(sourceVoxelId)
%                     sourceVoxelId = leafNode(this);
%                 end
                checkedBound = this.voxel_bounds(sourceVoxelId,:);
                queryVoxelBound = this.voxel_bounds(queryVoxelId,:);
                bool = all(bsxfun(@le, queryVoxelBound(:,1:3),checkedBound(:,4:6))&...
                           bsxfun(@ge, queryVoxelBound(:,4:6),checkedBound(:,1:3)),2);
                voxelNeighbourId = sourceVoxelId(bool);
            end
        end
        % Query equal adjacent voxels of the given voxel
        function voxelNeighbourId = queryEqualVoxelNeighbour(this, queryVoxelId, sourceVoxelId)
            if numel(queryVoxelId)>1
                warning('This function is used to query a adjacent neighbour of one voxel');
                return;
            else
                if isempty(sourceVoxelId)
                    sourceVoxelId = leafNode(this);
                end
                checkedBound = this.voxel_bounds(sourceVoxelId,:);
                queryVoxelBound = this.voxel_bounds(queryVoxelId,:);
                bool = all(bsxfun(@le,queryVoxelBound(:,1:2),checkedBound(:,4:5))&...
                           bsxfun(@ge,queryVoxelBound(:,4:5),checkedBound(:,1:2)),2)&...
                       all(bsxfun(@eq,queryVoxelBound(:,3),checkedBound(:,3))&...
                           bsxfun(@eq,queryVoxelBound(:,6),checkedBound(:,6)),2);
                voxelNeighbourId = sourceVoxelId(bool);
                voxelNeighbourId = setdiff(voxelNeighbourId, queryVoxelId);
            end
        end
        % Query lower adjacent voxels of the given voxel
        function voxelNeighbourId = queryLowerVoxelNeighbour(this, queryVoxelId, sourceVoxelId)
            if numel(queryVoxelId)>1
                warning('This function is used to query a adjacent neighbour of one voxel');
                return;
            else
                if isempty(sourceVoxelId)
                    sourceVoxelId = leafNode(this);
                end
                checkedBound = this.voxel_bounds(sourceVoxelId,:);
                queryVoxelBound = this.voxel_bounds(queryVoxelId,:);
                bool = all(bsxfun(@le, queryVoxelBound(:,1:2),checkedBound(:,4:5))&...
                           bsxfun(@ge, queryVoxelBound(:,4:5),checkedBound(:,1:2)),2)&...
                       all(bsxfun(@eq, queryVoxelBound(:,3),checkedBound(:,6)),2);
                voxelNeighbourId = sourceVoxelId(bool);
                voxelNeighbourId = setdiff(voxelNeighbourId, queryVoxelId);
            end
        end
        % Query equal and lower adjacent voxels of the given voxel
        function voxelNeighbourId = queryEqualLowerVoxelNeighbour(this, queryVoxelId, sourceVoxelId)
            if numel(queryVoxelId)>1
                warning('This function is used to query a adjacent neighbour of one voxel');
                return;
            else
                if isempty(sourceVoxelId)
                    sourceVoxelId = leafNode(this);
                end
                checkedBound = this.voxel_bounds(sourceVoxelId,:);
                queryVoxelBound = this.voxel_bounds(queryVoxelId,:);
                bool = all(bsxfun(@le,queryVoxelBound(:,1:2),checkedBound(:,4:5))&...
                           bsxfun(@ge,queryVoxelBound(:,4:5),checkedBound(:,1:2)),2)&...
                       all(bsxfun(@le,queryVoxelBound(:,3),checkedBound(:,6))&...
                           bsxfun(@gt,queryVoxelBound(:,6),checkedBound(:,3)),2);
                voxelNeighbourId = sourceVoxelId(bool);
                voxelNeighbourId = setdiff(voxelNeighbourId, queryVoxelId);
            end
        end
        % Query equal voxels of the given voxel
        function voxelNeighbourId = queryEqualVoxel(this, queryVoxelId, sourceVoxelId)
            if numel(queryVoxelId)>1
                warning('This function is used to query a adjacent neighbour of one voxel');
                return;
            else
                if isempty(sourceVoxelId)
                    sourceVoxelId = leafNode(this);
                end
                checkedBound = this.voxel_bounds(sourceVoxelId,:);
                queryVoxelBound = this.voxel_bounds(queryVoxelId,:);
                bool = all(bsxfun(@eq,queryVoxelBound(:,3),checkedBound(:,3))&...
                           bsxfun(@eq,queryVoxelBound(:,6),checkedBound(:,6)),2);
                voxelNeighbourId = sourceVoxelId(bool);
            end
        end
        % Window search
        function neighbour_voxel_ids = windowQueryVoxelNeighbour(this, query_voxel_ids, target_voxel_ids, windowSize)
            % Check input variables
            if isempty(query_voxel_ids)
                warning('The query voxel must not be empty');
                return;
            end
            if isempty(target_voxel_ids)
                warning('The query voxel must not be empty');
                return;
            end
                
            % Searching
            flag = false(numel(target_voxel_ids),1);
            target_voxel_bounds = this.voxel_bounds(target_voxel_ids,:);
            for i = 1:numel(query_voxel_ids)
                query_voxel_bounds = this.voxel_bounds(query_voxel_ids(i),:);
                
                mask = all(bsxfun(@le, query_voxel_bounds(:,1:3) - windowSize/2, target_voxel_bounds(:,4:6))&...
                       bsxfun(@ge, query_voxel_bounds(:,4:6) + windowSize/2, target_voxel_bounds(:,1:3)),2);
                neighbour_voxel_ids = target_voxel_ids(mask);
                neighbour_voxel_ids = setdiff(neighbour_voxel_ids, query_voxel_ids(i));
                
                mask = ismember(target_voxel_ids, neighbour_voxel_ids);
                flag(mask) = true;
            end
            neighbour_voxel_ids = target_voxel_ids(flag);      
        end
        
        % Box search
        function voxelNeighbourId = boxVoxelNeighbour(this, boxCoords, sourceVoxelId)
            if isempty(boxCoords) || size(boxCoords,2) < 6
                warning('The box must be defined');
                return;
            else
                checkedBound = this.voxel_bounds(sourceVoxelId,:);
                bool = all(bsxfun(@le, boxCoords(:,1:3), checkedBound(:,4:6))&...
                           bsxfun(@ge, boxCoords(:,4:6), checkedBound(:,1:3)),2);
                voxelNeighbourId = sourceVoxelId(bool);
            end
        end
    end
    % Computing feature of the voxel
    methods 
        % The normal OR the largest tangent vector OR the smoothest OR
        % residual  of the voxell
        function voxelFeature = voxelPCA(this, threshold, sourceVoxelId, varargin)
            %Feature = 1: normal vector
            %Feature = 2: a largest tangent vector
            %Feature = 3: The smoothest
            %Feature = 4: The residual
            %Feature = 5: All
            if isempty(varargin)
                warning('No plotting since not enough input variables');
                return;
            end
            IP = inputParser;
%             IP.addParamValue('numMinPts',3);
            IP.addParameter('feature',1);
            IP.parse(varargin{:})
            calVoxelFeature = IP.Results;
            % Check input
%             validateattributes(calVoxelFeature.numMinPts,{'numeric'},{'real','finite','nonnegative','nonnan'})
            if (calVoxelFeature.feature<1)||(calVoxelFeature.feature>5)
                calVoxelFeature.feature = 1; %Cal normal vector
            end
            %
            ref_vector = threshold.nz;
            if isempty(sourceVoxelId)
                voxel_Id = leafNode(this);
            else
                voxel_Id = sourceVoxelId;
            end
            voxelFeature = inf(numel(voxel_Id),7); %[id, normal/tangent,smooth/
            voxelFeature(:,1) = voxel_Id;
            for i=1:numel(voxel_Id)
                voxelPtsMark = this.voxel_ptc_ids == voxel_Id(i);
                voxelPts = this.Points(voxelPtsMark,:);
                if size(voxelPts,1) >= threshold.minNumPtsSurf
                    if calVoxelFeature.feature == 1 %Normal vector
                         normalvector = eigenspace(voxelPts,1);
                         if sign(dot(normalvector,ref_vector)) > 0
                             voxelFeature(i,2:4) = normalvector;
                         else
                             voxelFeature(i,2:4) = -normalvector;
                         end
                    elseif calVoxelFeature.feature == 2 %The tangent vector
                        voxelFeature(i,2:4) = eigenspace(voxelPts,2);
                    elseif calVoxelFeature.feature == 3 %The smoothest
                        voxelFeature(i,5) = eigenspace(voxelPts,3);
                    else % calCellFeature.feature == 4: Residual
                        [eigVector, eigValue, addResult]  = eigenspace(voxelPts,4);
                        %result = eigenspace(voxelPts,4);
                        voxelFeature(i,2:4) = eigVector(1,:);
                        voxelFeature(i,6)= addResult;
                        voxelFeature(i,7)= eigValue(2)/eigValue(3);
                    end
                end
            end
%             voxelFeature = voxelFeature(~isinf(voxelFeature(:,2)),:);
            if (calVoxelFeature.feature == 1)||(calVoxelFeature.feature == 2)
                voxelFeature(:,5:7) = [];
            elseif calVoxelFeature.feature == 3
                voxelFeature(:,6:7) = [];
                voxelFeature(:,2:4) = [];
            elseif calVoxelFeature.feature == 4
                voxelFeature(:,7) = [];
                voxelFeature(:,2:5) = [];
            else %feature  = 5
                voxelFeature(:,5) = [];
            end  
        end
        % Calculate a elevation gradient of the leaf voxel based on the points within the voxel 
        function voxelZGradient = elevGradient(this, varargin)
            %Feature = 1: Z difference
            %Feature = 2: Z gradient
            %Fearture = 3: slope
            %Feature = 4: Z difference and gradient
            IP = inputParser;
            IP.addParameter('voxelId',[]);
            IP.addParameter('feature',1);
            IP.parse(varargin{:})
            calVoxelZFeature = IP.Results;
            if (calVoxelZFeature.feature < 0)||(calVoxelZFeature.feature >4)
                calVoxelZFeature.feature = 1;
            end
            %
            if ~isempty(calVoxelZFeature.voxelId)
                sourceVoxelId = calVoxelZFeature.voxelId;
            else
                sourceVoxelId = leafNode(this);
            end
            voxelZGradient = zeros(numel(sourceVoxelId),4);
            voxelZGradient(:,1) = sourceVoxelId;
            for i=1:numel(sourceVoxelId)
                voxelPtsMark = this.voxel_ptc_ids == sourceVoxelId(i);
                voxelPts = this.Points(voxelPtsMark,:);
                if size(voxelPts,1) < 3
                    voxelZGradient(i,2) = inf;
                else
                    if calVoxelZFeature.feature == 1
                        voxelZGradient(i,2) = diff([min(voxelPts(:,3)),max(voxelPts(:,3))]);   
                    elseif calVoxelZFeature.feature == 2
                        [~,id_max] = max(voxelPts(:,3));
                        [~,id_min] = min(voxelPts(:,3));
%                         delta = voxelPts(id_max,:)-voxelPts(id_min,:);
%                         voxelZGradient(i,2) =sqrt((delta(3)/delta(1))^2+(delta(3)/delta(2))^2);                       
                        voxelZGradient(i,2) = (voxelPts(id_max,3)-voxelPts(id_min,3))/...
                                              norm(voxelPts(id_max,1:2)-voxelPts(id_min,1:2));
                    elseif calVoxelZFeature.feature == 3%Slope
                        ref_vector = [0,0,1];
                        normalvector = eigenspace(voxelPts,1);
                        if sign(dot(normalvector,ref_vector)) < 0
                            normalvector = -normalvector;
                        end
                        voxelZGradient(i,2) = tan(dot(normalvector,ref_vector));                   
                    else
                        [~,id_max] = max(voxelPts(:,3));
                        [~,id_min] = min(voxelPts(:,3));
                        voxelZGradient(i,2) = diff([voxelPts(id_min,3),voxelPts(id_max,3)]);
                        voxelZGradient(i,3) = diff([voxelPts(id_min,3),voxelPts(id_max,3)])/...
                                              norm(voxelPts(id_max,1:2)-voxelPts(id_min,1:2));     
                        ref_vector = [0,0,1];
                        normalvector = eigenspace(voxelPts,1);
                        if sign(dot(normalvector,ref_vector)) < 0
                            normalvector = -normalvector;
                        end
                        voxelZGradient(i,4) = tan(dot(normalvector,ref_vector));                  
                                          
                    end
                end
            end
            voxelZGradient = voxelZGradient(~isinf(voxelZGradient(:,2)),:);
            if calVoxelZFeature.feature ~= 4
                voxelZGradient(:,3:4) = [];
            end
        end
        %Compute a normal vector, tangent vector and deepest slope
        function result = voxelGeometry(this, varargin)
            IP = inputParser;
            IP.addParameter('voxelId',[]);
            IP.addParameter('minNumPts',1);
            IP.parse(varargin{:})
            calVoxelGeometry = IP.Results;
            if isempty(calVoxelGeometry.voxelId)
                voxelId = leafNode(this);
            else
                voxelId = calVoxelGeometry.voxelId;
            end
            if isempty(calVoxelGeometry.minNumPts)||(calVoxelGeometry.minNumPts<3)
                minNumPoint = 3;
            else
                minNumPoint = calVoxelGeometry.minNumPts;
            end
            nz = [0,0,1];
            % Computer
            result = zeros(size(voxelId,1),8);
            result(:,1) = voxelId;
            for i=1:size(voxelId,1)
                voxelPtsMark = this.voxel_ptc_ids == voxelId(i);
                voxelPts = this.Points(voxelPtsMark,:);
                if size(voxelPts,1)<minNumPoint
                    result(i,2:end) = inf;
                else
                    normalvector = eigenspace(voxelPts,1);
                    cosAngle = dot(nz,normalvector);
                    if sign(cosAngle)<0
                        normalvector = -normalvector;
                    end
                    t = normalvector - nz/cosAngle;
                    t = t/norm(t);
                    %Update normal and tangent vectors
                    result(i,2:4) = normalvector;
                    result(i,5:7) = t;
                    %Compute the deepest slope                
                    P = mean(voxelPts,1);
                    [~,idx] = max(voxelPts(:,3));
                    P1 = voxelPts(idx,:);
                    [~,idx] = min(voxelPts(:,3));
                    P2 = voxelPts(idx,:);
                    % Point projection
%                     P11 = P+dot(t,P1-P)*t;
%                     P22 = P+dot(t,P2-P)*t;
                    % Plane intersection
                    delta1 = -(dot(nz,P)-P1(3))/dot(nz,t);
                    P11 = P1+delta1*t;
                    delta2 = -(dot(nz,P)-P2(3))/dot(nz,t);
                    P22 = P2+delta2*t;
                    result(i,8) = (P11(3)-P22(3))/norm(P11(1:2)-P22(1:2));
                end
            end
        end
    end % End method
    methods %Visualization
        function plotAll(this, varargin)
            if isempty(varargin)
                warning('No plotting since not enough input variables');
                return;
            end
            IP = inputParser;
            IP.addParameter('octreeDepth',[]);
            IP.addParameter('voxelId',[]);
            IP.addParameter('typePlotVoxel',2);
            IP.addParameter('typePlotPoint',0);
            IP.addParameter('color',[]);
            IP.addParameter('fill',1); % 1: filled cubic; 0; wireframe
            IP.parse(varargin{:})
            plotOption = IP.Results;
            %
            if isempty(plotOption.octreeDepth)&&isempty(plotOption.voxelId)
                voxelId = leafNode(this);
            elseif ~isempty(plotOption.octreeDepth)&& isempty(plotOption.voxelId)
                if (plotOption.octreeDepth<0)||(plotOption.octreeDepth>max(this.depth))
                    plotOption.octreeDepth = max(this.depth);
                end
                voxelId = this.depth == plotOption.octreeDepth;
            else
                voxelId = plotOption.voxelId;
            end
            %
            if (plotOption.typePlotVoxel<0)||(plotOption.typePlotVoxel>2)
                typePlotVoxel = 2; %Both empty and full voxel
            else
                typePlotVoxel = plotOption.typePlotVoxel; %Either empty "0" or full "1"
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
            if typePlotVoxel~=2 %Either empty "0" or full "1"
                voxelId = voxelId(this.voxel_prop(voxelId)==typePlotVoxel);
            end
            listVoxel = (1:this.voxel_count)'; 
            voxelId = listVoxel(voxelId);
            % Plot
            plotVoxel(this, voxelId, typePlotVoxel, color, fill)
            if (typePlotVoxel~=0)&&(typePlotPoint==1)&& (fill==0)
                plotPoint(this,voxelId, color) 
            end
        end
        % Plot voxel   
        function plotVoxel(this, voxelId, typePlotVoxel, color, fill)
            hold on;
            if isempty(color)
               color = lines(numel(voxelId));
            else
                color = repmat(color,[numel(voxelId),1]);
            end
            for i = 1:numel(voxelId)
                voxelMinMax = this.voxel_bounds(voxelId(i),:);
                if fill==1 % Plot filled voxel   
                    face1 = cat(1,voxelMinMax([1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3]));
                    face2 = cat(1,voxelMinMax([1 2 6; 1 5 6; 4 5 6; 4 2 6; 1 2 6]));
                    face3 = cat(1,voxelMinMax([1 2 3; 1 2 6; 4 2 6; 4 2 3; 1 2 3]));
                    face4 = cat(1,voxelMinMax([4 2 3; 4 2 6; 4 5 6; 4 5 3; 4 2 3]));
                    face5 = cat(1,voxelMinMax([4 5 3; 4 5 6; 1 5 6; 1 5 3; 4 5 3]));
                    face6 = cat(1,voxelMinMax([1 5 3; 1 5 6; 1 2 6; 1 2 3; 1 5 3]));
                    if typePlotVoxel==2
                        if this.voxel_prop(voxelId(i))==1
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
                else %No fill voxel
                    voxelEdge = cat(1, voxelMinMax([...
                        1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3;...
                        1 2 6; 4 2 6; 4 5 6; 1 5 6; 1 2 6; 1 2 3]),...
                        nan(1,3), voxelMinMax([4 2 3; 4 2 6]),...
                        nan(1,3), voxelMinMax([4 5 3; 4 5 6]),...
                        nan(1,3), voxelMinMax([1 5 3; 1 5 6]));
                    if typePlotVoxel==2
                        if this.voxel_prop(voxelId(i))==1
                           plot3(voxelEdge(:,1),voxelEdge(:,2),voxelEdge(:,3),...
                               'Color',color(i,:),'LineWidth', 2);
                        else
                            plot3(voxelEdge(:,1),voxelEdge(:,2),voxelEdge(:,3),...
                               'Color','k','LineWidth', 2);
                        end
                    else
                        plot3(voxelEdge(:,1),voxelEdge(:,2),voxelEdge(:,3),...
                               'Color',color(i,:),'LineWidth', 2);
                    end
                end
            end
        end 
        % Plot points
        function plotPoint(this,voxelId,color)
            hold on
            doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
            if isempty(color)
               color = lines(numel(voxelId));
            else
                color = repmat(color,[numel(voxelId),1]);
            end
            for i = 1:numel(voxelId)
                this.Points(this.voxel_ptc_ids(voxelId(i)).id,:)
%                 doplot3(this.Points(this.voxel_ptc_ids==voxelId(i),:),'.','Color',color(i,:),'MarkerSize',5)
                doplot3(this.Points(this.voxel_ptc_ids(voxelId(i)).id,:),'.','Color',color(i,:),'MarkerSize',5)
            end 
        end
        %Plot voxels according to its depth and order in each layer
        function plotOrderVoxel(this)
            maxDepth = max(this.depth);
            numVoxel = this.voxel_count;
            voxelId = 1:numVoxel;
            h = zeros(numVoxel,1);
            count = 0;
            hold all;
            for i=1:maxDepth+1
                bool = this.depth==(i-1);
                voxelId_level = voxelId(bool);
                colors = lines(numel(voxelId_level));
                for j=1:numel(voxelId_level)
                    count = count+1;
                    voxelMinMax = this.voxel_bounds(voxelId_level(j),:);
                    voxelEdge = cat(1, voxelMinMax([...
                                    1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3;...
                                    1 2 6; 4 2 6; 4 5 6; 1 5 6; 1 2 6; 1 2 3]),...
                                    nan(1,3), voxelMinMax([4 2 3; 4 2 6]),...
                                    nan(1,3), voxelMinMax([4 5 3; 4 5 6]),...
                                    nan(1,3), voxelMinMax([1 5 3; 1 5 6]));
                    h(count) = plot3(voxelEdge(:,1),voxelEdge(:,2),voxelEdge(:,3));
                    set(h(count),'Color',colors(j,:),'LineWidth', 2)
                    pause
                end
                
            end
        end
        % Plot voxel coloring by its attribute
        function plotVoxelAttribute(this, varargin)
            IP = inputParser;
            IP.addParameter('voxelId',[]);
            IP.addParameter('attribute',[]);
            IP.addParameter('colormap','hsv');
            IP.addParameter('fill',1); % 1: filled cubic; 0; wireframe
            IP.parse(varargin{:})
            plotOption = IP.Results;
            if isempty(plotOption.voxelId)||isempty(plotOption.attribute)
                return;
            end
            %
            voxelId = plotOption.voxelId;
            attribute = plotOption.attribute;
            cmap = colormap(plotOption.colormap);
            normalizeAtt = (attribute-min(attribute))/(max(attribute)-min(attribute));
            rangeAtt = linspace(min(normalizeAtt),max(normalizeAtt),size(cmap,1)+1);
            for i=1:numel(rangeAtt)-1
                bool = (normalizeAtt>=rangeAtt(i))&(normalizeAtt<rangeAtt(i+1));
                curVoxelId = voxelId(bool);
                if ~isempty(curVoxelId)
                    for j=1:numel(curVoxelId)
                        voxelMinMax = this.voxel_bounds(curVoxelId(j),:);
                        if plotOption.fill==1
                            face1 = cat(1,voxelMinMax([1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3]));
                            face2 = cat(1,voxelMinMax([1 2 6; 1 5 6; 4 5 6; 4 2 6; 1 2 6]));
                            face3 = cat(1,voxelMinMax([1 2 3; 1 2 6; 4 2 6; 4 2 3; 1 2 3]));
                            face4 = cat(1,voxelMinMax([4 2 3; 4 2 6; 4 5 6; 4 5 3; 4 2 3]));
                            face5 = cat(1,voxelMinMax([4 5 3; 4 5 6; 1 5 6; 1 5 3; 4 5 3]));
                            face6 = cat(1,voxelMinMax([1 5 3; 1 5 6; 1 2 6; 1 2 3; 1 5 3]));
                            %
                            fill3(face1(:,1),face1(:,2),face1(:,3),cmap(i,:))
                            fill3(face2(:,1),face2(:,2),face2(:,3),cmap(i,:))
                            fill3(face3(:,1),face3(:,2),face3(:,3),cmap(i,:))
                            fill3(face4(:,1),face4(:,2),face4(:,3),cmap(i,:))
                            fill3(face5(:,1),face5(:,2),face5(:,3),cmap(i,:))
                            fill3(face6(:,1),face6(:,2),face6(:,3),cmap(i,:))
                        else
                            voxelEdge = cat(1, voxelMinMax([...
                                1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3;...
                                1 2 6; 4 2 6; 4 5 6; 1 5 6; 1 2 6; 1 2 3]),...
                                nan(1,3), voxelMinMax([4 2 3; 4 2 6]),...
                                nan(1,3), voxelMinMax([4 5 3; 4 5 6]),...
                                nan(1,3), voxelMinMax([1 5 3; 1 5 6]));
                            plot3(voxelEdge(:,1),voxelEdge(:,2),voxelEdge(:,3),...
                               'Color',cmap(i,:),'LineWidth', 2);
                        end
                    end         
                end
            end
%             h = colorbar;
            cbh = colorbar('YGrid','on');
            num_space = size(get(cbh,'ytick'),2)-1;
            colorLabel = cell(1,num_space+1);
            plot_dist = round(linspace(min(attribute),max(attribute),num_space+1)*100)/100;
            for j=1:numel(plot_dist)
                colorLabel{j} = sprintf('%0.2f',plot_dist(j));
            end
            set(cbh,'yticklabel',colorLabel,'FontSize',10,'FontName','Time New Roman');
        end
    end %End a method visualization
end