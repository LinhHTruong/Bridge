function [w_ptc_center, w_ptc_normal] = wPCA(varargin)

IP = inputParser;
%The Data set
IP.addParameter('ptc',[]);
% Define initial normal vector
IP.addParameter('ini_normal_vector', [0.5, 0.5, 0.5]/norm([0.5, 0.5, 0.5]));
% Define angle threshold
IP.addParameter('angle_threshold',5*pi/180);
% maximum iteraction
IP.addParameter('max_iters',100);
% Define a minimum distance threshold
IP.addParameter('min_pts', 5);  

IP.parse(varargin{:})
% Assign the value
ptc = IP.Results.ptc;
ini_normal_vector = IP.Results.ini_normal_vector;
cosin_threshold = cos(IP.Results.angle_threshold);
max_iters = IP.Results.max_iters;
min_pts = IP.Results.min_pts;

% Normalize the initial vector
ini_normal_vector = ini_normal_vector./norm(ini_normal_vector);
% Compute the normal vector
if size(ptc, 1) > min_pts
    % Initial center and normal
    ptc_center = mean(ptc, 1);
    if isempty(ini_normal_vector)
        ptc_normal = eigenspace(ptc, 1);
    else
        ptc_normal = ini_normal_vector;
    end
    
    % Iterations
    flag = true;
    iters = 0;
    while flag && (iters <= max_iters)

        % Establish weight for each point based on a distance
        dist_ptc_surface = abs(dist_3Dpoints_3Dplane(ptc, [ptc_center, ptc_normal]));
        if all(dist_ptc_surface)
            w_ptc_center = ptc_center;
            w_ptc_normal = ptc_normal;
            return
        end
        ptc_weight = exp(-dist_ptc_surface.^2/(mean(dist_ptc_surface)/5)^2);
        
        % Computer new center and normal
        w_ptc = bsxfun(@times, ptc, ptc_weight);
        w_ptc_center = sum(w_ptc, 1)./sum(ptc_weight);
        w_ptc_cov = (ptc - w_ptc_center)'*(ptc - w_ptc_center)/sum(ptc_weight);
        
        [eig_vects, eig_vals] = eig(w_ptc_cov);
        [~, mask]= sort(diag(eig_vals));
        w_ptc_normal = eig_vects(:, mask(1))';

        % Check deviation of the normal
        if abs(dot(w_ptc_normal, ptc_normal)) <  cosin_threshold % Angle deviation is still large
            ptc_center = w_ptc_center;
            ptc_normal = w_ptc_normal;
            iters = iters + 1;
        else
            flag = false;
        end
    end
    
else
    w_ptc_center = [inf, inf, inf];
    w_ptc_normal = [inf, inf, inf];
end
