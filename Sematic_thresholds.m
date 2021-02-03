classdef Sematic_thresholds < handle
    properties
        sampling_step
        nz
        cell_size
        min_num_cell
        min_num_pts
        max_angle
        max_residual
        max_distance
        range_search_radius
        min_dist_2_planes
        hor_plane_max_angle
        voxel_size
   
    end
    methods
        function this = Sematic_thresholds(varargin)           
            IP = inputParser;
            IP.addParameter('sampling_step',0.005);
            IP.addParameter('nz',[0,0,1]);
            % The cell size is used in subdivide the space
            IP.addParameter('cell_size',1.0); 
            % Minimum number of the cells in the region
            IP.addParameter('min_num_cell',10);
            % The minimum number of points within the cell or neighbour for computing eigenvector 
            IP.addParameter('min_num_pts',5);
            % Define angle between the cell or region
            IP.addParameter('max_angle',15);  
            % A residual to consider the cell as the smooth surface represent to perfect surfaces of the 
            % structure, which can be a seeding cell
            IP.addParameter('max_residual', 0.01);
            % A distance to consider for the voxel region growing
            IP.addParameter('max_distance', 0.01); 
            % Range search to find neighbour
            IP.addParameter('range_search_radius', 2.5);
            % The smalest size of the structural components, which imply the distance between 
            % different two parallel cell is no smaller than 0.2m
            IP.addParameter('min_dist_2_planes', 0.20);
            % The maximum angle between the plane and oz, where the plane
            % can considered non vertical plane
            IP.addParameter('hor_plane_max_angle', 60);
            % The voxel size is used in subdivide the space
            IP.addParameter('voxel_size',0.1);
            
            IP.parse(varargin{:})
            
            % Assign the value
            this.nz = IP.Results.nz;
            this.sampling_step = IP.Results.sampling_step;
            this.cell_size = IP.Results.cell_size;
            this.min_num_cell = IP.Results.min_num_cell;
            this.min_num_pts = IP.Results.min_num_pts;
            this.max_angle = IP.Results.max_angle;
            this.max_residual = IP.Results.max_residual;
            this.max_distance = IP.Results.max_distance;
            this.range_search_radius = IP.Results.range_search_radius;
            this.min_dist_2_planes = IP.Results.min_dist_2_planes;
            this.hor_plane_max_angle = IP.Results.hor_plane_max_angle;
            this.voxel_size = IP.Results.voxel_size;
        end
    end
end