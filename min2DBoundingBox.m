function min_bb = min2DBoundingBox(X)
% compute the minimum bounding box of a set of 2D points
%   Use:   boundingBox = min2DBoundingBox(point_matrix)
%
% Input:  Nx2 matrix containing the [x,y] coordinates of n points
%         *** there must be at least 3 points which are not collinear
% output: 4x2 matrix containing the coordinates of the bounding box corners
%
% Example : generate a random set of point in a randomly rotated rectangle
%     n = 50000;
%     t = pi*rand(1);
%     X = [cos(t) -sin(t) ; sin(t) cos(t)]*[7 0; 0 2]*rand(2,n);
%     X = [X  20*(rand(2,1)-0.5)];  % add an outlier
% 
%     tic
%     c = minBoundingBox(X);
%     toc
% 
%     figure(42);
%     hold off,  plot(X(1,:),X(2,:),'.')
%     hold on,   plot(c(1,[1:end 1]),c(2,[1:end 1]),'r')
%     axis equal
% compute the convex hull (CH is a 2*k matrix subset of X)
% compute the convex hull (CH is a 2*k matrix subset of X)
k = convhull(X(1,:),X(2,:));
CH = X(:,k);
% compute the angle to test, which are the angle of the CH edges as:
%   "one side of the bounding box contains an edge of the convex hull"
E = diff(CH,1,2);           % CH edges
T = atan2(E(2,:),E(1,:));   % angle of CH edges (used for rotation)
T = unique(mod(T,pi/2));    % reduced to the unique set of first quadrant angles
% create rotation matrix which contains
% the 2x2 rotation matrices for *all* angles in T
% R is a 2n*2 matrix
R = cos( reshape(repmat(T,2,2),2*length(T),2) ... % duplicate angles in T
       + repmat([0 -pi ; pi 0]/2,length(T),1));   % shift angle to convert sine in cosine
% rotate CH by all angles
RCH = R*CH;
% compute border size  [w1;h1;w2;h2;....;wn;hn]
% and area of bounding box for all possible edges
bsize = max(RCH,[],2) - min(RCH,[],2);
area  = prod(reshape(bsize,2,length(bsize)/2));
% find minimal area, thus the index of the angle in T 
[~,i] = min(area);
% compute the bound (min and max) on the rotated frame
Rf    = R(2*i+[-1 0],:);   % rotated frame
bound = Rf * CH;           % project CH on the rotated frame
bmin  = min(bound,[],2);
bmax  = max(bound,[],2);
% compute the corner of the bounding box
Rf = Rf';
bb(:,4) = bmax(1)*Rf(:,1) + bmin(2)*Rf(:,2);
bb(:,1) = bmin(1)*Rf(:,1) + bmin(2)*Rf(:,2);
bb(:,2) = bmin(1)*Rf(:,1) + bmax(2)*Rf(:,2);
bb(:,3) = bmax(1)*Rf(:,1) + bmax(2)*Rf(:,2);
bb = bb';
% compute features of the bounding box
bb_edge = bb(2:3,:) - bb(1:2,:);
bb_edge_length = sqrt(sum(bsxfun(@power, bb_edge, 2),2));
[~, long_edge_id] = max(bb_edge_length);
[~, short_edge_id] = min(bb_edge_length);

min_bb = struct('vertices',[],'area',[],'long_edge',[],'long_edge_vector',[],'short_edge',[],'short_edge_vector',[]);
min_bb.vertices = bb;
min_bb.long_edge = bb_edge_length(long_edge_id);
min_bb.long_edge_vector = bb_edge(long_edge_id,:)/norm(bb_edge(long_edge_id,:));
min_bb.short_edge = bb_edge_length(short_edge_id);
min_bb.short_edge_vector = bb_edge(short_edge_id,:)/norm(bb_edge(short_edge_id,:));
min_bb.area = min_bb.long_edge*min_bb.short_edge;


