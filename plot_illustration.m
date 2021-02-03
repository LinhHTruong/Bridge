%% Plot the cell + plane extraction
threshold = THRESHOLD;
Tree = OQTR;
bandwidth = 0.5*0.125*threshold.min_dist_2_planes;
cell_id = 2264;%2266;%2276;
cell_ptc_ids = Tree.cell_pts(cell_id).id;
cell_ptc_xyz = Tree.pts(cell_ptc_ids,1:3);

% Adjust points
mask = cell_ptc_xyz(:,3) >= 1.0;
cell_ptc_xyz(mask,3) = cell_ptc_xyz(mask,3) - 0.6;


% Extract the points in the cell after removing coincide points in the
% vertical direction
[~,ia,~] = unique(cell_ptc_xyz(:,1:2), 'rows');
cell_filter_ptc_xyz = cell_ptc_xyz(ia,:);

% Kernel density estimation
cell_filter_ptc_range = [min(cell_filter_ptc_xyz(:,3)), max(cell_filter_ptc_xyz(:,3))];
no_kde_pt = max(100,ceil(diff(cell_filter_ptc_range)/bandwidth));
[fi,zi,~] = ksdensity(cell_filter_ptc_xyz(:,3),'npoints',no_kde_pt,'bandwidth',bandwidth,'Kernel','epanechnikov');

% Calculate the first derivative
st_fi = deriv(zi, fi);
st_fi = deriv(zi, st_fi);

% n =length(zi);
% st_fi(1) = (fi(2) - fi(1))/(zi(2) - zi(1));
% st_fi(n) = (fi(n) - fi(n-1))/(zi(n) - zi(n-1));
% for j = 2:n-1
%     st_fi(j)=(fi(j+1)-fi(j))./(zi(j+1) - zi(j));
% %         d(j)=(y(j+1)-y(j-1))./(x(j+1) - x(j-1));
% end
    
    
peak_shape = peak_shape_width(zi, fi);
mask = peak_shape(:,3) + peak_shape(:,2) > 0;
peak_shape = peak_shape(mask,:);
% Sort from heigh to low
[~, mask] = sort(peak_shape(:,1), 'descend');
peak_shape = peak_shape(mask,:);
clear cell_filter_ptc_xyz
peak_shape(:,4) = peak_shape(:,2) + peak_shape(:,3);
% peak_shape(:,[2,3]) = [];
    
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 900 900]);
set(gca,'DataAspectRatio',[1 1 1]);

colormap('lines');
cmap = colormap;

%
% plot3(mean(cell_ptc_xyz(:,1))-fi/10, ones(numel(zi),1)*mean(cell_ptc_xyz(:,2)), zi,'-', 'linewidth', 2,'color','r')
% Set data and kde
fi = 0.5*(fi/max(fi));
kde_x = mean(cell_ptc_xyz(:,1))-fi;
kde_y = ones(1,numel(zi))*mean(cell_ptc_xyz(:,2));
kde_z = zi;
kde = [kde_x',kde_y',kde_z'];

% first derivative kde
st_fi = 0.4*(st_fi/max(st_fi));
st_kde_x = mean(cell_ptc_xyz(:,1))-st_fi;
st_kde_y = ones(1,numel(zi))*mean(cell_ptc_xyz(:,2));
st_kde_z = zi;
st_kde = [st_kde_x',st_kde_y',st_kde_z'];

% Rotate
kde_cent = mean(kde,1);
kde = kde - kde_cent;

st_kde_cent = mean(st_kde,1);
st_kde = st_kde - st_kde_cent;

alpha = 30;
Rz = [cos(alpha)    -sin(alpha)     0;...
      sin(alpha)     cos(alpha)     0;...
      0              0              1];  
kde = kde*Rz;
kde = kde + kde_cent;

st_kde = st_kde*Rz;
st_kde = st_kde + kde_cent;

% Set label
axis_interval = [0.2, 0.2, 0.1];
x_range = [min(min(cell_ptc_xyz(:,1)), min([kde(:,1);st_kde(:,1)])), max(max(cell_ptc_xyz(:,1)), max([kde(:,1);st_kde(:,1)]))];
y_range = [min(min(cell_ptc_xyz(:,2)), min([kde(:,2);st_kde(:,2)])),max(max(cell_ptc_xyz(:,2)), max([kde(:,2);st_kde(:,2)]))];
z_range = [min(min(cell_ptc_xyz(:,3)), min([kde(:,3);st_kde(:,3)])), max(max(cell_ptc_xyz(:,3)), max([kde(:,3);st_kde(:,3)]))];

x_range = [floor(x_range(1)*10)/10, ceil(x_range(2)*10)/10];
y_range = [floor(y_range(1)*10)/10, ceil(y_range(2)*10)/10];
z_range = [floor(z_range(1)*10)/10, ceil(z_range(2)*10)/10];


xyz_cent = [mean(x_range), mean(y_range), mean(z_range)];


num_space = ([diff(x_range), diff(y_range), diff(z_range)]./axis_interval);
% num_space = ceil([diff(x_range), diff(y_range), diff(z_range)]./axis_interval);
% num_space = 2*ceil(num_space/2);

x_range = xyz_cent(1)-axis_interval(1)*num_space(1)/2:axis_interval(1):xyz_cent(1)+axis_interval(1)*num_space(1)/2;
y_range = xyz_cent(2)-axis_interval(2)*num_space(2)/2:axis_interval(2):xyz_cent(2)+axis_interval(2)*num_space(2)/2;
z_range = xyz_cent(3) -axis_interval(3)*num_space(3)/2:axis_interval(3):xyz_cent(3)+axis_interval(3)*num_space(3)/2;

% 

% view(105,18)
% Plot original data

subplot(1,2,1)
hold all

% Kde
h1 = plot3(kde(:,1), kde(:,2), kde(:,3),'-', 'linewidth', 2,'color','r');
% kde axis
plot3([max(kde(:,1));max(kde(:,1))], [max(kde(:,2));max(kde(:,2))],[min(kde(:,3)); max(kde(:,3))], '-', 'linewidth', 1,'color',[98,98,98]/255)

% kde second derivative

h2 = plot3(st_kde(:,1), st_kde(:,2), st_kde(:,3),'-', 'linewidth', 2,'color','b');

% Plot the line
% plot_style = {'b-', 'k-', 'r-', 'm-', 'c-'};
patch_ptc_ids = [];
for i = 1:size(peak_shape,1)
    for j = 2:3
        if j == 2
            z1_1 = peak_shape(i,1) - peak_shape(i,j);
        else
            z1_1 = peak_shape(i,1) + peak_shape(i,j);
        end
        % find fi
        mask = kde(:,3) == z1_1;
        kde_x1 = kde(mask,1);
        kde_y1 = kde(mask,2);
        mask = st_kde(:,3) == z1_1;
        st_kde_x1 = st_kde(mask,1);
        st_kde_y1 = st_kde(mask,2);
        plot3([kde_x1; st_kde_x1],[kde_y1;st_kde_y1],[z1_1, z1_1], '-', 'linewidth', 2,'color',cmap(i+2,:))
        
%         Plot the point on KDE
        plot3(kde_x1, kde_y1, z1_1, 'marker', 'd', 'markersize', 8, 'markerface', cmap(2,:))
        
        
    end
    
    % Plot the points of the local surface
%     z1 = peak_shape(i,1) - peak_shape(i,2);
%     z2 = peak_shape(i,1) + peak_shape(i,3);
%     mask = (z1 <= cell_ptc_xyz(:,3))&(cell_ptc_xyz(:,3) <= z2);
%     plot3(cell_ptc_xyz(mask,1), cell_ptc_xyz(mask,2), cell_ptc_xyz(mask,3),'.', 'markersize', 10,'color',cmap(i+2,:))

    patch_ptc_ids = [patch_ptc_ids;find(mask)];
end
% Plot original points
plot3(cell_ptc_xyz(:,1), cell_ptc_xyz(:,2), cell_ptc_xyz(:,3),'.', 'markersize', 5,'color',[128,128,128]/255)

% Plot unsurface points
% mask = ismember([1:size(cell_ptc_xyz,1)],patch_ptc_ids);
% plot3(cell_ptc_xyz(~mask,1), cell_ptc_xyz(~mask,2), cell_ptc_xyz(~mask,3),'.', 'markersize', 3,'color',[128,128,128]/255)



% lg = legend(h1,{sprintf('KDE')},'Location','north');
lg = legend([h1, h2],[{sprintf('KDE')}, {sprintf('Second derivative of KDE')}],'Location','northwest');
lg.Orientation = 'horizontal';
lg.FontSize = 18;
lg.FontName = 'Time New Roman';
lg.NumColumns = 1;
lg.Box = 'off';
    
    
% Set lims
ax = gca;
ax.TickLength = [0.01, 0.01];
ax.XAxis.TickLabelFormat = '%.1f';
ax.YAxis.TickLabelFormat = '%.1f';
ax.ZAxis.TickLabelFormat = '%.1f';
ax.FontName = 'Time New Roman';
ax.FontSize = 18;
ax.XLim = [min(x_range), max(x_range)];
ax.YLim = [min(y_range), max(y_range)];
ax.ZLim = [min(z_range), max(z_range)];
ax = gca;
ax.XTick = x_range;
ax.YTick = y_range;
ax.ZTick = z_range;

txt = ax.XLabel;
txt.String = 'x';
txt.FontSize = 18;

txt = ax.YLabel;
txt.String = 'y';
txt.FontSize = 18;

txt = ax.ZLabel;
txt.String = 'z';
txt.FontSize = 18;


grid on
view(-115,10)

% Save file

% file_name = sprintf('Origin_ptc_cell_%d',cell_id);
% file_name = sprintf('Origin_ptc_cell_%d_KDE',cell_id);
file_name = sprintf('Origin_ptc_cell_%d_KDE_ndKDE',cell_id);
% file_name = sprintf('Origin_ptc_cell_%d_patches',cell_id);
saveas(gcf,strcat(file_name,'.png'))
F = getframe(gcf);
imwrite(F.cdata,strcat(file_name,'.tif'),'Resolution',300) 

%% Plot road curb extraction progress
% Use the version 4, for Germany bridge, cell size = 0.5
curb_cell_id = 2264;
% Road
st_cell_surf;
st_cell_pts_ids;
st_surf_cent_1 = st_cell_surf(1:3);
st_surf_cent_2 = st_surf_cent_1 + 0.15*st_cell_surf(4:6);

%footpath
nd_cell_surf
nd_cell_pts_ids;
nd_surf_cent_1 = nd_cell_surf(1:3);
nd_surf_cent_2 = nd_surf_cent_1 + 0.15*nd_cell_surf(4:6);

% Cell points
cell_pts_ids = Tree.cell_pts(curb_cell_id).id;
cell_pts_xyz = Tree.pts(cell_pts_ids,1:3);

% Unassign points
road_curb_cell_unassign_pts_ids = setdiff(cell_pts_ids, union(st_cell_pts_ids, nd_cell_pts_ids));
                

% Extract the points between the surfaces
[ptc_local_ids, ~] = extract_ptc_2_planes(cell_pts_xyz, st_cell_surf, nd_cell_surf);
road_curb_cell_pts_ids = cell_pts_ids(ptc_local_ids);

% Remove points within two surface
road_curb_cell_remain_pts_ids = setdiff(road_curb_cell_pts_ids, union(st_cell_pts_ids, nd_cell_pts_ids));
                

% Filter the points on the road curb
mask = ismember(cell_pts_ids, road_curb_cell_remain_pts_ids);
road_curb_cell_remain_pts_xyz = cell_pts_xyz(mask,:);
pts_local_id = road_curb_ransac_filter(road_curb_cell_remain_pts_xyz, threshold);

% Final road curb points
road_curb_cell_final_pts_ids = road_curb_cell_remain_pts_ids(pts_local_id);
road_curb_cell_final_pts_xyz = road_curb_cell_remain_pts_xyz(pts_local_id,1:3);
                

%%    
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 900 900]);
set(gca,'DataAspectRatio',[1 1 1]);
hold all
colormap('lines');
cmap = colormap;

% Get the points within the cell
mask = cell_pts_xyz(:,3) >= st_cell_surf(3) - 0.01;
cell_sub_pts_xyz = cell_pts_xyz(mask,:);

% Set label
axis_interval = [0.1, 0.1, 0.05];
x_range = [min(cell_sub_pts_xyz(:,1)), max(cell_sub_pts_xyz(:,1))];
y_range = [min(cell_sub_pts_xyz(:,2)), max(cell_sub_pts_xyz(:,2))];
z_range = [min(cell_sub_pts_xyz(:,3)), max(max(cell_sub_pts_xyz(:,3)), nd_surf_cent_2(3))];

x_range = [floor(x_range(1)*10)/10, ceil(x_range(2)*10)/10];
y_range = [floor(y_range(1)*10)/10, ceil(y_range(2)*10)/10];
z_range = [floor(z_range(1)*10)/10, ceil(z_range(2)*10)/10];

xyz_cent = [mean(x_range), mean(y_range), mean(z_range)];
num_space = ([diff(x_range), diff(y_range), diff(z_range)]./axis_interval);
% num_space = ceil([diff(x_range), diff(y_range), diff(z_range)]./axis_interval);
% num_space = 2*ceil(num_space/2);

x_range = xyz_cent(1)-axis_interval(1)*num_space(1)/2:axis_interval(1):xyz_cent(1)+axis_interval(1)*num_space(1)/2;
y_range = xyz_cent(2)-axis_interval(2)*num_space(2)/2:axis_interval(2):xyz_cent(2)+axis_interval(2)*num_space(2)/2;
z_range = xyz_cent(3) -axis_interval(3)*num_space(3)/2:axis_interval(3):xyz_cent(3)+axis_interval(3)*num_space(3)/2;

% Plot points
% h1 = plot3(cell_sub_pts_xyz(:,1), cell_sub_pts_xyz(:,2), cell_sub_pts_xyz(:,3),'.', 'markersize', 9,'color',[128,128,128]/255)

% Plot unassign points 
% plot3(Tree.pts(road_curb_cell_unassign_pts_ids,1), Tree.pts(road_curb_cell_unassign_pts_ids,2), Tree.pts(road_curb_cell_unassign_pts_ids,3),'.', 'markersize', 9,'color',cmap(1,:));
% h2 = plot3(nan, nan, nan,'.', 'markersize', 30,'color',cmap(1,:));


% Plot the road points
plot3(Tree.pts(st_cell_pts_ids,1), Tree.pts(st_cell_pts_ids,2), Tree.pts(st_cell_pts_ids,3),'.', 'markersize', 9,'color',cmap(2,:));
h3 = plot3(nan, nan, nan,'.', 'markersize', 30,'color',cmap(2,:));

% % Plot surface
% plot3(st_surf_cent_1(1), st_surf_cent_1(2), st_surf_cent_1(3),'.', 'markersize', 30,'color',cmap(1,:))
% mArrow3(st_surf_cent_1,st_surf_cent_2, 'stemWidth',0.0025,'facealpha',1, 'color', cmap(1,:));
% bound = [min(Tree.pts(st_cell_pts_ids,:)), max(Tree.pts(st_cell_pts_ids,:))];
% corners = bound([1, 2, 3; 4, 2, 3; 4, 5, 6;1, 5, 3; 1, 2, 3]);
% corners(:,3) = st_cell_surf(3);
% face1 = corners;
% h4 = fill3(face1(:,1),face1(:,2),face1(:,3),cmap(1,:));
% set(h4,'facealpha',.2)
% h4.EdgeColor = cmap(1,:);

% Plot the footpath points
plot3(Tree.pts(nd_cell_pts_ids,1), Tree.pts(nd_cell_pts_ids,2), Tree.pts(nd_cell_pts_ids,3),'.', 'markersize', 9,'color',cmap(5,:));
h5 = plot3(nan, nan, nan,'.', 'markersize', 30,'color',cmap(5,:));

% 
% % Plot surface
% h8 = plot3(nd_surf_cent_1(1), nd_surf_cent_1(2), nd_surf_cent_1(3),'.', 'markersize', 30,'color',cmap(3,:));
% mArrow3(nd_surf_cent_1,nd_surf_cent_2, 'stemWidth',0.0025,'facealpha',1, 'color', cmap(3,:));
% bound = [min(Tree.pts(nd_cell_pts_ids,:)), max(Tree.pts(nd_cell_pts_ids,:))];
% corners = bound([1, 2, 3; 4, 2, 3; 4, 5, 6;1, 5, 3; 1, 2, 3]);
% corners(:,3) = nd_cell_surf(3);
% face1 = corners;
% h9 = fill3(face1(:,1),face1(:,2),face1(:,3),cmap(3,:));
% set(h9,'facealpha',.2)
% h9.EdgeColor = cmap(3,:);

% Plot remaining points
% plot3(road_curb_cell_remain_pts_xyz(:,1),road_curb_cell_remain_pts_xyz(:,2), road_curb_cell_remain_pts_xyz(:,3),'.', 'markersize', 10,'color',[128,128,128]/255);
% h10 = plot3(nan, nan, nan,'.', 'markersize', 30,'color',[128,128,128]/255);

% 
% % Plot final points
plot3(road_curb_cell_final_pts_xyz(:,1), road_curb_cell_final_pts_xyz(:,2), road_curb_cell_final_pts_xyz(:,3),'.', 'markersize', 10,'color',cmap(6,:));
h11 = plot3(nan, nan, nan,'.', 'markersize', 30,'color',cmap(6,:));
% 

% Set lims
ax = gca;
ax.TickLength = [0.01, 0.01];
ax.XAxis.TickLabelFormat = '%.2f';
ax.YAxis.TickLabelFormat = '%.2f';
ax.ZAxis.TickLabelFormat = '%.2f';
ax.FontName = 'Time New Roman';
ax.FontSize = 18;
ax.XLim = [min(x_range), max(x_range)];
ax.YLim = [min(y_range), max(y_range)];
ax.ZLim = [min(z_range), max(max(z_range),nd_surf_cent_2(3))];
ax = gca;
ax.XTick = x_range;
ax.YTick = y_range;
ax.ZTick = z_range;

txt = ax.XLabel;
txt.String = 'x';
txt.FontSize = 24;

txt = ax.YLabel;
txt.String = 'y';
txt.FontSize = 24;

txt = ax.ZLabel;
txt.String = 'z';
txt.FontSize = 24;
grid on
view(-73,24)

% Legend
% lg = legend([h3, h5, h2],[{sprintf('Road points')}, {sprintf('Footpath points')}, {sprintf('Unassigned points')}]);
% lg.ItemTokenSize(1) = 10;
% lg = legend([h4, h9],[{sprintf('Road')}, {sprintf('Footpath')}]);
% lg.ItemTokenSize(1) = 60;
% lg = legend([h3, h5, h10],[{sprintf('Road points')}, {sprintf('Footpath points')}, {sprintf('Candidate points')}]);
% lg.ItemTokenSize(1) = 10;
lg = legend([h3, h5, h11],[{sprintf('Road points')}, {sprintf('Footpath points')}, {sprintf('Road curb points')}]);
lg.ItemTokenSize(1) = 10;

% [lg, lg_obj] = legend([h2, h3, h5],[{sprintf('Unassigned points')}, {sprintf('Road points')}, {sprintf('Footpath points')}]);
% lg_obj1 = findobj(lg_obj, 'type', 'line'); 
% % lg_obj1.Markersize = 20;
% lg_obj1(2).MarkerSize = 30
% lg_obj1(4).MarkerSize = 20
% lg_obj1(6).MarkerSize = 20

lg.Location = 'north';
lg.Orientation = 'horizontal';
lg.FontSize = 24;
lg.FontName = 'Time New Roman';
lg.NumColumns = 3;
lg.Box = 'off';
% 
% 
% lg.Position = lg_position
% file_name = sprintf('Road_curb_cell_%d_points',curb_cell_id);
% file_name = sprintf('Road_curb_cell_%d_Fitting_Surface',curb_cell_id);
% file_name = sprintf('Road_curb_cell_%d_Vertical points',curb_cell_id);
file_name = sprintf('Road_curb_cell_%d_Final Road Curb',curb_cell_id);
F = getframe(gcf);
imwrite(F.cdata,strcat(file_name,'.tif'),'Resolution',300)
saveas(gcf,strcat(file_name,'.png'))
 
%% Road parapet

% Version 4 with cell_size = 1.0/Sub_Conrete bridge


% Footpath points
i = 1;
footpath_id = footpath_ids(i);
footpath_cell_ids = vertcat(supstruct(1).Component(footpath_id).cell.ids);
footpath_cell_pts_ids = vertcat(supstruct(1).Component(footpath_id).cell.ptc_ids);
footpath_cell_pts_xyz = Tree.pts(footpath_cell_pts_ids,1:3);
footpath_pts_range = [min(footpath_cell_pts_xyz); max(footpath_cell_pts_xyz)];
% Road points
mask = cellfun(@(s) contains(s, 'Road Surface'), supstruct(1).Description);
road_id = find(mask);
road_cell_ids = vertcat(supstruct(1).Component(road_id).cell.ids);
road_cell_pts_ids = vertcat(supstruct(1).Component(road_id).cell.ptc_ids);
road_cell_pts_xyz = Tree.pts(road_cell_pts_ids,1:3);
road_pts_range = [min(road_cell_pts_xyz); max(road_cell_pts_xyz)];

% Road points
mask = cellfun(@(s) contains(s, 'Road Curb 1'), supstruct(1).Description);
road_curb_id = find(mask);
road_curb_cell_ids = vertcat(supstruct(1).Component(road_curb_id).cell.ids);
road_curb_cell_pts_ids = vertcat(supstruct(1).Component(road_curb_id).cell.ptc_ids);
road_curb_cell_pts_xyz = Tree.pts(road_curb_cell_pts_ids,1:3);

% rot_road_mbb = min2DBoundingBox(road_cell_pts_xyz(:,1:2)');
% rot_road_mbb.vertices(:,3) = mean(road_cell_pts_xyz(:,3));
% rot_road_mbb.long_edge_vector(:,3) = mean(road_cell_pts_xyz(:,3));
% rot_road_mbb.short_edge_vector(:,3) = mean(road_cell_pts_xyz(:,3));
% rot_road_mbb.long_edge_vector = rot_road_mbb.long_edge_vector/norm(rot_road_mbb.long_edge_vector);
% rot_road_mbb.short_edge_vector = rot_road_mbb.short_edge_vector/norm(rot_road_mbb.short_edge_vector);


road_cent = mean(road_surface.mbb.vertices);
road_dir = road_surface.mbb.long_edge_vector;


assign_ptc_ids = [footpath_cell_pts_ids;road_cell_pts_ids;road_curb_cell_pts_ids];

% Get color
ids = randperm(256);
ids = ids(1:64);

%% Update road because the algorithm calculates incorrect for this case
% rot_road_mbb = min2DBoundingBox(road_cell_pts_xyz(:,1:2)');
% rot_road_mbb.vertices(:,3) = mean(road_cell_pts_xyz(:,3));
% rot_road_mbb.long_edge_vector(:,3) = 0;
% rot_road_mbb.short_edge_vector(:,3) = 0;
% rot_road_mbb.long_edge_vector = rot_road_mbb.short_edge_vector/norm(rot_road_mbb.short_edge_vector);
% rot_road_mbb.short_edge_vector = rot_road_mbb.long_edge_vector/norm(rot_road_mbb.long_edge_vector);
% road_surface.mbb = rot_road_mbb; % use to call distance in main codes

%%
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 900 900]);
set(gca,'DataAspectRatio',[1 1 1]);
hold all
colormap('lines');
cmap = colormap;

% Plot the points within the cell
% 
% for j = 1:numel(parapet_cell_ids)
%     parapet_cell_id = parapet_cell_ids(j);
%     parapet_cell_pts_ids = Tree.cell_pts(parapet_cell_id).id;
%     parapet_cell_pts_ids = setdiff(parapet_cell_pts_ids, assign_ptc_ids);
%     parapet_cell_pts_xyz = Tree.pts(parapet_cell_pts_ids,1:3);
%     
%     plot3(parapet_cell_pts_xyz(:,1), parapet_cell_pts_xyz(:,2), parapet_cell_pts_xyz(:,3),'.', 'markersize', 5,'color',cmap(ids(j),:));
%  
% end


% % Plot the the parapet

% for j = 1:length(parapet.cell)
% %     % Points of the cells
% %     parapet_cell_pts_ids = parapet.cell(j).ptc_ids;
% %     parapet_cell_pts_xyz = Tree.pts(parapet_cell_pts_ids,1:3);
% %     
% %     % Plot points in cells
% %     plot3(parapet_cell_pts_xyz(:,1), parapet_cell_pts_xyz(:,2), parapet_cell_pts_xyz(:,3),'.', 'markersize', 5,'color',cmap(ids(j),:));
%     
%     % Plot cell centers
%     parapet_cell_bounds_cent = parapet.cell(j).surface_features;
%     plot3(parapet_cell_bounds_cent(:,1), parapet_cell_bounds_cent(:,2), max(footpath_pts_range(:,3))+0.1,'.', 'markersize', 30,'color',cmap(ids(j),:));
%  
% end
% 

% Plot final parapet
for count = 1:length(final_parapet)
    
    parapet_ptc_ids = vertcat(final_parapet(count).cell.ptc_ids);
    parapet_ptc_xyz = Tree.pts(parapet_ptc_ids,1:3);
    plot3(parapet_ptc_xyz(:,1), parapet_ptc_xyz(:,2), parapet_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(count + 5,:));
%     
end
            
            
            
% Plot footpath
plot3(footpath_cell_pts_xyz(:,1), footpath_cell_pts_xyz(:,2), footpath_cell_pts_xyz(:,3),'.', 'markersize', 20,'color',[128,128,128]/255);


% Plot road curb
plot3(road_curb_cell_pts_xyz(:,1), road_curb_cell_pts_xyz(:,2), road_curb_cell_pts_xyz(:,3),'.', 'markersize', 20,'color',[168,168,168]/255);


% Plot road
plot3(road_cell_pts_xyz(:,1), road_cell_pts_xyz(:,2), road_cell_pts_xyz(:,3),'.', 'markersize', 20,'color',[198,198,198]/255);
% road_cent(3) = road_cent(3) + 0.1;
% plot3(road_cent(1), road_cent(2), road_cent(3),'.', 'markersize', 30,'color',cmap(1,:))
% road_cent_2 = road_cent + 3*road_dir;
% road_cent_2(3) = road_cent(3);
% mArrow3(road_cent, road_cent_2, 'stemWidth',0.03,'facealpha',1, 'color', cmap(7,:));

axis off
view(21,33)

% file_name = sprintf('Road_parapet_points');
% file_name = sprintf('Road_parapet_candiate_points');
% file_name = sprintf('Road_parapet_Cell_Projection');
file_name = sprintf('Road_parapet_Final');
F = getframe(gcf);
imwrite(F.cdata,strcat(file_name,'.tif'),'Resolution',300)
saveas(gcf,strcat(file_name,'.png'))
%% Plot the distance distribution
no_kde_pt = 100;
[fi,zi,~] = ksdensity(dist_parapet_road,'npoints',no_kde_pt,'bandwidth',0.5*bin_width_dist,'Kernel','epanechnikov');

% Calculate the first derivative
st_fi = deriv(zi, fi);
st_fi = deriv(zi, st_fi);

peak_shape = peak_shape_width(zi, fi);
mask = peak_shape(:,3) + peak_shape(:,2) > 0;
peak_shape = peak_shape(mask,:);
% Sort from heigh to low
[~, mask] = sort(peak_shape(:,1), 'descend');
peak_shape = peak_shape(mask,:);
clear cell_filter_ptc_xyz
peak_shape(:,4) = peak_shape(:,2) + peak_shape(:,3);
% peak_shape(:,[2,3]) = [];
    
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 900 900]);
set(gca,'DataAspectRatio',[1 1 1]);
hold all
colormap('lines');
cmap = colormap;

% Plot cell_center
cell_bounds_proj_cent = parapet_cell_bounds_proj_cent;
cell_bounds_proj_cent(:,3) = max(footpath_pts_range(:,3))+0.05;

% Plot footpath
plot3(footpath_cell_pts_xyz(:,1), footpath_cell_pts_xyz(:,2), footpath_cell_pts_xyz(:,3),'.', 'markersize', 20,'color',[128,128,128]/255);
% 
% 
% Plot road curb
plot3(road_curb_cell_pts_xyz(:,1), road_curb_cell_pts_xyz(:,2), road_curb_cell_pts_xyz(:,3),'.', 'markersize', 20,'color',[168,168,168]/255);
% 
% 
% % Plot road
plot3(road_cell_pts_xyz(:,1), road_cell_pts_xyz(:,2), road_cell_pts_xyz(:,3),'.', 'markersize', 20,'color',[198,198,198]/255);
% road_cent(3) = road_cent(3) + 0.1;
% plot3(road_cent(1), road_cent(2), road_cent(3),'.', 'markersize', 30,'color',cmap(1,:))
% road_cent_2 = road_cent + 3*road_dir;
% road_cent_2(3) = road_cent(3);
% mArrow3(road_cent, road_cent_2, 'stemWidth',0.03,'facealpha',1, 'color', cmap(7,:));
road_cent(3) = max(road_cell_pts_xyz(:,3))+0.05
road_y_min = min(road_cell_pts_xyz(:,2));
road_y_max = max(road_cell_pts_xyz(:,2));
y_pos = linspace(road_y_min, road_y_max, 4);
for j = 1:numel(y_pos)-1
    
    pos_1 = [road_cent(1), y_pos(j), road_cent(3)];
    pos_2 = pos_1 + (y_pos(j+1) - y_pos(j))*road_dir;
    % Start points
    plot3(pos_1(1), pos_1(2), pos_1(3),'.', 'markersize', 30,'color',cmap(1,:))
    
    % Direction
    
    mArrow3(pos_1, pos_2, 'stemWidth',0.03,'facealpha',1, 'color', cmap(7,:));


end

for j = 1:length(parapet.cell)
    % Plot cell centers
    plot3(cell_bounds_proj_cent(j,1), cell_bounds_proj_cent(j,2), cell_bounds_proj_cent(j,3),'.', 'markersize', 30,'color',cmap(ids(j),:));
 
end


% Set data and kde
fi = 3*(fi/max(fi));
kde_x =  road_surface.cent(1) - zi + 0;
kde_y = ones(1,numel(zi))*road_surface.cent(2);
kde_z = mean(road_cent(:,3))+fi;
kde = [kde_x',kde_y',kde_z'];

% first derivative kde
st_fi = 3*(st_fi/max(st_fi));
st_kde_x = road_surface.cent(1) - zi;
st_kde_y = ones(1,numel(zi))*road_surface.cent(2);
st_kde_z = mean(road_cent(:,3))+st_fi;
st_kde = [st_kde_x',st_kde_y',st_kde_z'];

% Rotate
kde_cent = mean(kde,1);
kde = kde - kde_cent;

st_kde_cent = mean(st_kde,1);
st_kde = st_kde - st_kde_cent;

alpha = 0;
Rz = [cos(alpha)    -sin(alpha)     0;...
      sin(alpha)     cos(alpha)     0;...
      0              0              1];  
kde = kde*Rz;
kde = kde + kde_cent;

st_kde = st_kde*Rz;
st_kde = st_kde + kde_cent;

% Kde
h1 = plot3(kde(:,1), kde(:,2), kde(:,3),'-', 'linewidth', 2,'color','r');
% kde axis
plot3([min(kde(:,1));max(kde(:,1))], [min(kde(:,2));min(kde(:,2))],[min(kde(:,3)); min(kde(:,3))], '-', 'linewidth', 1,'color',[98,98,98]/255)

% kde second derivative

% h2 = plot3(st_kde(:,1), st_kde(:,2), st_kde(:,3),'-', 'linewidth', 2,'color','b');


axis off
view(21,33)
file_name = sprintf('Road_parapet_KDE_multi_cell');
F = getframe(gcf);
imwrite(F.cdata,strcat(file_name,'.tif'),'Resolution',300)
saveas(gcf,strcat(file_name,'.png'))


%% Plot cell_based region growing
% Back and forward : Subdata set, cell_size = 1
% cell_ids = [63, 65, 71, 62, 64, 70, 78, 80, 86, 79, 81, 87];
% Back and forward : Subdata set_ro15, cell_size = 1
cell_ids = [44,45,48,49,74, 75,78,79,76,77,80,81, 108, 130, 132]; 
    
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 900 900]);
set(gca,'DataAspectRatio',[1 1 1]);
hold all
colormap('lines');
colormap('hsv');
cmap = colormap;

color_ids = 1:4:12*numel(cell_ids);
ptc_range = [];
% Plot the points of the cell
all_cell_ptc_xyz = [];
for i = 1:numel(cell_ids)
    % Get points within the cell
    cell_id = cell_ids(i);
    mask = region_info(:,1) == cell_id;
    cell_peak_ids = region_info(mask,[1,2]);
    mask = ismember(cells_plane.cell_ids(:,[1,2]), cell_peak_ids, 'rows');
    cell_ptc_ids = cells_plane.peak_info(mask).ptc_ids;
    cell_ptc_xyz = Tree.pts(cell_ptc_ids,1:3);
    ptc_range = [ptc_range ; min(cell_ptc_xyz, [], 1),max(cell_ptc_xyz, [], 1)];
    all_cell_ptc_xyz = [all_cell_ptc_xyz;cell_ptc_xyz]; 
    
    % Plot points
%     plot3(cell_ptc_xyz(:,1), cell_ptc_xyz(:,2), cell_ptc_xyz(:,3),'.', 'markersize', 4,'color',cmap(color_ids(i),:));
%     plot3(cell_ptc_xyz(:,1), cell_ptc_xyz(:,2), cell_ptc_xyz(:,3),'.', 'markersize', 10,'color',[128 128, 128]/255);
end
% % % axis off

% Plot points redering according to elevation
num_interval = 16;
color_ids = ceil(linspace(1, 255, num_interval));
elev_pos = linspace(min(all_cell_ptc_xyz(:,3)), max(all_cell_ptc_xyz(:,3)), num_interval+1);

for i = 1:num_interval
    mask = (elev_pos(i) <= all_cell_ptc_xyz(:,3))&(all_cell_ptc_xyz(:,3) <= elev_pos(i+1));
    plot3(all_cell_ptc_xyz(mask,1), all_cell_ptc_xyz(mask,2), all_cell_ptc_xyz(mask,3),'.', 'markersize', 4,'color',cmap(color_ids(i),:));
end

%     
% Plot cells according to the region
% mask = ismember(region_info(:,1), cell_ids);
% region_cell_ids = region_info(mask,:);
% region_ids = unique(region_cell_ids(:,3));
% for i = 1:numel(region_ids)
%     % Get points within the cell
%     mask = region_cell_ids(:,3) == region_ids(i);
%     cell_peak_ids = region_cell_ids(mask,[1,2]);
%     mask = ismember(cells_plane.cell_ids(:,[1,2]), cell_peak_ids, 'rows');
%     cell_ptc_ids = vertcat(cells_plane.peak_info(mask).ptc_ids);
%     cell_ptc_xyz = Tree.pts(cell_ptc_ids,1:3);
%     ptc_range = [ptc_range ; min(cell_ptc_xyz, [], 1),max(cell_ptc_xyz, [], 1)];
%     
%     % Plot points
%     plot3(cell_ptc_xyz(:,1), cell_ptc_xyz(:,2), cell_ptc_xyz(:,3),'.', 'markersize', 10,'color',cmap(color_ids(i),:));
% end

% Plot results of back/forward
% for i = 1:numel(region_ids)
% %     % Get cell_peak
% %     mask = region_cell_ids(:,3) == region_ids(i);
% %     cell_peak_ids = region_cell_ids(mask,[1,2]);
% %     
%     % Get points
%     cell_peak = vertcat(cells_region(region_ids(i)).cell.id);
%     mask = ismember(cell_peak(:,1), cell_ids);
%     cell_ptc_ids = vertcat(cells_region(region_ids(i)).cell(mask).ptc_ids);
%     cell_ptc_xyz = Tree.pts(cell_ptc_ids,1:3);
%     ptc_range = [ptc_range ; min(cell_ptc_xyz, [], 1),max(cell_ptc_xyz, [], 1)];
%     
%     % Plot
%     plot3(cell_ptc_xyz(:,1), cell_ptc_xyz(:,2), cell_ptc_xyz(:,3),'.', 'markersize', 10,'color',cmap(color_ids(i),:));
% end
   
    
% Set range
ptc_range = [min(ptc_range(:,1:3),[],1), max(ptc_range(:,4:6),[], 1)];
% Set label
axis_interval = [0.5, 0.5, 0.25];

ptc_range(1:3) = floor(ptc_range(1:3)*10)/10;
ptc_range(4:6) = ceil(ptc_range(4:6)*10)/10;

xyz_cent = 0.5*(ptc_range(1:3) + ptc_range(4:6));

num_space = ceil([ptc_range(4:6) - ptc_range(1:3)]./axis_interval);
% num_space = ceil([diff(x_range), diff(y_range), diff(z_range)]./axis_interval);
% num_space = 2*ceil(num_space/2);

x_range = xyz_cent(1)-axis_interval(1)*num_space(1)/2:axis_interval(1):xyz_cent(1)+axis_interval(1)*num_space(1)/2;
y_range = xyz_cent(2)-axis_interval(2)*num_space(2)/2:axis_interval(2):xyz_cent(2)+axis_interval(2)*num_space(2)/2;
z_range = xyz_cent(3) -axis_interval(3)*num_space(3)/2:axis_interval(3):xyz_cent(3)+axis_interval(3)*num_space(3)/2;


% Set lims
ax = gca;
ax.TickLength = [0.01, 0.01];
ax.XAxis.TickLabelFormat = '%.2f';
ax.YAxis.TickLabelFormat = '%.2f';
ax.ZAxis.TickLabelFormat = '%.2f';
ax.FontName = 'Time New Roman';
ax.FontSize = 24;
ax.XLim = [min(x_range), max(x_range)];
ax.YLim = [min(y_range), max(y_range)];
ax.ZLim = [min(z_range), max(z_range)];
ax = gca;
ax.XTick = x_range;
ax.YTick = y_range;
ax.ZTick = z_range;

txt = ax.XLabel;
txt.String = 'x';
txt.FontSize = 24;

txt = ax.YLabel;
txt.String = 'y';
txt.FontSize = 24;

txt = ax.ZLabel;
txt.String = 'z';
txt.FontSize = 24;
grid on

% Legend
% lg = legend([h3, h5, h2],[{sprintf('Road points')}, {sprintf('Footpath points')}, {sprintf('Unassigned points')}]);
% lg.ItemTokenSize(1) = 10;
% lg = legend([h4, h9],[{sprintf('Road')}, {sprintf('Footpath')}]);
% lg.ItemTokenSize(1) = 60;
% lg = legend([h3, h5, h10],[{sprintf('Road points')}, {sprintf('Footpath points')}, {sprintf('Candidate points')}]);
% lg.ItemTokenSize(1) = 10;
% lg = legend([h3, h5, h11],[{sprintf('Road points')}, {sprintf('Footpath points')}, {sprintf('Road curb points')}]);
% lg.ItemTokenSize(1) = 10;

% [lg, lg_obj] = legend([h2, h3, h5],[{sprintf('Unassigned points')}, {sprintf('Road points')}, {sprintf('Footpath points')}]);
% lg_obj1 = findobj(lg_obj, 'type', 'line'); 
% % lg_obj1.Markersize = 20;
% lg_obj1(2).MarkerSize = 30
% lg_obj1(4).MarkerSize = 20
% lg_obj1(6).MarkerSize = 20

% lg.Location = 'north';
% lg.Orientation = 'horizontal';
% lg.FontSize = 24;
% lg.FontName = 'Time New Roman';
% lg.NumColumns = 3;
% lg.Box = 'off';





    
view(35,35)
% 
% 
% lg.Position = lg_position
file_name = sprintf('Back_Forward_Points_Elevation');
% file_name = sprintf('Back_Forward_Cell_points');
% file_name = sprintf('Back_Forward_Segment_Points');
% file_name = sprintf('Back_Forward_Final_Points');
% file_name = sprintf('Back_Forward_Final_Points_Boundary_Filter');
F = getframe(gcf);
imwrite(F.cdata,strcat(file_name,'.tif'),'Resolution',300)
saveas(gcf,strcat(file_name,'.png'))



%% Merge region
source_level_id = 5;
level_id = 6;
count = length(Cells_Region(source_level_id).cell);
num_cells = length(Cells_Region(level_id).cell);
for j = 1:num_cells
    Cells_Region(source_level_id).cell(count + 1).id = Cells_Region(level_id).cell(j).id;
    Cells_Region(source_level_id).cell(count + 1).ptc_ids = Cells_Region(level_id).cell(j).ptc_ids;
    Cells_Region(source_level_id).cell(count + 1).surface_features = Cells_Region(level_id).cell(j).surface_features;

    count = count + 1;
end
    


    
%% Plot section - top and bottom surface
ptc = PTC.xyz;
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 900 900]);
set(gca,'DataAspectRatio',[1 1 1]);
hold all
% colormap('lines');
colormap('jet');
cmap = colormap;



% Plot points redering according to elevation
num_interval = 32;
color_ids = ceil(linspace(1, 64, num_interval));
elev_pos = linspace(min(ptc(:,3)), max(ptc(:,3)), num_interval+1);
% 
% for i = 1:num_interval
%     mask = (elev_pos(i) <= ptc(:,3))&(ptc(:,3) <= elev_pos(i+1));
%     if (1.52<=elev_pos(i))&(elev_pos(i+1)<=1.74)
%         plot3(ptc(mask,1), ptc(mask,2), ptc(mask,3),'.', 'markersize', 20,'color',cmap(color_ids(i),:));
%     else
%         plot3(ptc(mask,1), ptc(mask,2), ptc(mask,3),'.', 'markersize', 3,'color',cmap(color_ids(i),:));
%     end
% end
% Create the colar bar
% colormap(cmap(color_ids,:));
% colorLabel = cell(1,num_interval+1);
% % for j=1:2:num_color+1
% 
% for j = 1:num_interval + 1
%     if mod(j,4) == 1 
%         colorLabel{j} = sprintf('%.2f',elev_pos(j));
%     else
%         colorLabel{j} = '';
%     end
% end
% cbh = colorbar;
% cbh.Position = [0.875,0.275, 0.04, 0.4];
% cbh.Label.String = 'Elevation (m)';
% cbh.FontSize = 20;
% cbh.Ticks = linspace(0,1,num_interval+1);
% cbh.TickLabels = colorLabel;
% cbh.FontName = 'Time New Roman';
% cbh.TickLength = 0;
% % cbh.Label.String = 'Deformation (m)';
% cbh.Label.Position = [cbh.Position(1)+0.5*cbh.Position(3),1.15,0];
% cbh.Label.Rotation = 0;
% cbh.Label.HorizontalAlignment = 'center';
% cbh.Label.VerticalAlignment = 'top';

% Plot all points
% plot3(ptc(:,1), ptc(:,2), ptc(:,3),'.', 'markersize', 10,'color',[168,168,168]/255);

% Plot the road surface
assign_ptc_ids = [];
road_cell_ptc_ids = vertcat(SuperStructure.Component(1).cell.ptc_ids);
road_cell_ptc_xyz = OQTR.pts(road_cell_ptc_ids,:);
plot3(road_cell_ptc_xyz(:,1), road_cell_ptc_xyz(:,2), road_cell_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(5,:));
h1 = plot3(nan, nan, nan,'.', 'markersize', 30,'color',cmap(5,:));

assign_ptc_ids= [assign_ptc_ids; road_cell_ptc_ids];
% Plot footpath
fpath_cell_ptc_ids = vertcat(SuperStructure.Component(2).cell.ptc_ids);
fpath_cell_ptc_xyz = OQTR.pts(fpath_cell_ptc_ids,:);
plot3(fpath_cell_ptc_xyz(:,1), fpath_cell_ptc_xyz(:,2), fpath_cell_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(20,:));
h2 = plot3(nan, nan, nan,'.', 'markersize', 30,'color',cmap(20,:));
assign_ptc_ids= [fpath_cell_ptc_ids; road_cell_ptc_ids];

% Plot footpath
fpath_cell_ptc_ids = vertcat(SuperStructure.Component(3).cell.ptc_ids);
fpath_cell_ptc_xyz = OQTR.pts(fpath_cell_ptc_ids,:);
plot3(fpath_cell_ptc_xyz(:,1), fpath_cell_ptc_xyz(:,2), fpath_cell_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(35,:));
h3 = plot3(nan, nan, nan,'.', 'markersize', 30,'color',cmap(35,:));

assign_ptc_ids= [fpath_cell_ptc_ids; road_cell_ptc_ids];


% Plot road direction
road_dir = Road_Surface.mbb.long_edge_vector;
road_cent = mean(Road_Surface.mbb.vertices);
road_cent(3) = max(road_cell_ptc_xyz(:,3))+0.05;
road_x_min = min(road_cell_ptc_xyz(:,1));
road_x_max = max(road_cell_ptc_xyz(:,1));
x_pos = linspace(road_x_min, road_x_max, 4);
for j = 1:numel(x_pos)-1
    
    pos_1 = [x_pos(j), road_cent(2), road_cent(3)];
    pos_2 = pos_1 + (x_pos(j+1) - x_pos(j))*road_dir;
    % Start points
    plot3(pos_1(1), pos_1(2), pos_1(3),'.', 'markersize', 30,'color',cmap(50,:))
    
    % Direction
    mArrow3(pos_1, pos_2, 'stemWidth',0.04,'facealpha',1, 'color', cmap(55,:));
end

axis off
view(-70,25)

% file_name = sprintf('Points of Section');
% file_name = sprintf('Top surfaces');
file_name = sprintf('Classification');
% file_name = sprintf('Back_Forward_Final_Points');
% file_name = sprintf('Back_Forward_Final_Points_Boundary_Filter');
F = getframe(gcf);
imwrite(F.cdata,strcat(file_name,'.tif'),'Resolution',300)
saveas(gcf,strcat(file_name,'.png'))

%% Plot bottom surfaces
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 900 900]);
set(gca,'DataAspectRatio',[1 1 1]);
hold all
% colormap('lines');
colormap('jet');
cmap = colormap;


% Plot the region
num_region = length(Cells_Region);
regions_pts = [];
for region_id = 1:num_region
    cell_ptc_ids = vertcat(Cells_Region(region_id).cell.ptc_ids);
    cell_ptc_xyz = OQTR.pts(cell_ptc_ids,1:3);
    % Plot region
    plot3(cell_ptc_xyz(:,1), cell_ptc_xyz(:,2), cell_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(region_id*12,:));
%     num_cells = length(Cells_Region(region_id).cell);
%     for j = 1:num_cells
%         % For each cell
%         cell_ptc_ids = Cells_Region(region_id).cell(j).ptc_ids;
%         cell_ptc_xyz = OQTR.pts(cell_ptc_ids,1:3);
%         
%         % Plot each cells
%         color_id = j - floor(j/size(cmap,1))*size(cmap,1)+1;
%         plot3(cell_ptc_xyz(:,1), cell_ptc_xyz(:,2), cell_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(color_id,:));
%     end
    
end

% Plot road surface
% colormap('hsv');
% cmap = colormap;
% num_cells = length(SuperStructure.Component(1).cell);
% for j = 1:num_cells
%     % For each cell
%     cell_ptc_ids = SuperStructure.Component(1).cell(j).ptc_ids;
%     cell_ptc_xyz = OQTR.pts(cell_ptc_ids,1:3);
% 
%     % Plot each cells
%     color_id = j - floor(j/size(cmap,1))*size(cmap,1)+1;
%     plot3(cell_ptc_xyz(:,1), cell_ptc_xyz(:,2), cell_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(color_id,:));
% end
road_cell_ptc_ids = vertcat(SuperStructure.Component(1).cell.ptc_ids);
road_cell_ptc_xyz = OQTR.pts(road_cell_ptc_ids,:);
plot3(road_cell_ptc_xyz(:,1), road_cell_ptc_xyz(:,2), road_cell_ptc_xyz(:,3),'.', 'markersize', 20,'color',[198,198,198]/255);





axis off
view(-70,25)
% file_name = sprintf('Bottom surfaces');
file_name = sprintf('Bottom_surfaces_road');
% file_name = sprintf('Back_Forward_Final_Points');
% file_name = sprintf('Back_Forward_Final_Points_Boundary_Filter');
F = getframe(gcf);
imwrite(F.cdata,strcat(file_name,'.tif'),'Resolution',300)
saveas(gcf,strcat(file_name,'.png'))


%% Plot identify the bottom surface
bin_width = 0.5*threshold.min_dist_2_planes;
num_bin = ceil((max(dist_region_road_surface(:,1))-min(dist_region_road_surface(:,1)))/bin_width);
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 900 900]);
% set(gca,'DataAspectRatio',[1 1 1]);
hold all
% colormap('lines');
colormap('jet');
cmap = colormap;

% Plot hist based distance
% hist(dist_region_road_surface(:,1), num_bin);
[n, xout] = hist(dist_region_road_surface(:,1), num_bin);
bin_width = diff(xout);

color = cmap(ceil(linspace(1, 63, numel(n))),:);

for i = 1:numel(n)
    h1 = bar(xout(i), n(i), 'barwidth', bin_width(1) );
    if i == 1, hold on, end
    h1.FaceColor = color(i,:);
    


end
% h1 = bar(xout, n)

% Add KDE
no_kde_pt = 100;
[fi,zi,~] = ksdensity(dist_region_road_surface(:,1),'npoints',no_kde_pt,'bandwidth',0.5*bin_width(1),'Kernel','epanechnikov');

plot(zi, 100*fi, '-', 'linewidth', 2,'color','b');



% Set label
% x1 = min(min(xout) - bin_width(1), min(zi));
% x2 = max(max(xout) + - bin_width(1), max(zi));
% x_num_space = ceil((x2-x1)/(0.5*bin_width(1)));
% 
% x_range = x1:bin_width(1)/2:x1+(x_num_space)*bin_width(1)/2;
x_range = [min(zi), xout, max(zi)];
y_interval = 50;
y_num_space = ceil(max(n)/y_interval);
y_range = [0 : y_interval: y_num_space*y_interval];




% Set lims
ax = gca;
ax.TickLength = [0.01, 0.01];
ax.XAxis.TickLabelFormat = '%.2f';
ax.YAxis.TickLabelFormat = '%.0f';
ax.FontName = 'Time New Roman';
ax.FontSize = 24;
ax.XLim = [min(x_range), max(x_range)];
ax.YLim = [min(y_range), max(y_range)];
ax.XTick = x_range;
ax.YTick = y_range;

txt = ax.XLabel;
txt.String = 'Distances from cells to a road surface (m)';
txt.FontSize = 24;

txt = ax.YLabel;
txt.String = 'Number of cells';
txt.FontSize = 24;

grid on

% color = cmap(1:25:64,:)
% h1 = findobj(gca,'Type','patch');
% h1.FaceColor = color;
% h1.EdgeColor = 'w';
% h1.LineWidth = 5;
file_name = sprintf('Distance_bottom_surface_Road');
F = getframe(gcf);
imwrite(F.cdata,strcat(file_name,'.tif'),'Resolution',300)
saveas(gcf,strcat(file_name,'.png'))

%% Plot distribution of the region id in the largest peaks
[fi,zi,~] = ksdensity(dist_region_road_surface(:,1),'npoints',100,'bandwidth',0.5*threshold.min_dist_2_planes,'Kernel','epanechnikov');
peak_shape = peak_shape_width(zi, fi);
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 900 900]);
% set(gca,'DataAspectRatio',[1 1 1]);
hold all
% colormap('lines');
colormap('jet');
cmap = colormap;

% Plot hist based distance
count = 2;
mask = (peak_shape(count,1) - peak_shape(count,2) <= dist_region_road_surface(:,1))&...
           (dist_region_road_surface(:,1) <= peak_shape(count,1) + peak_shape(count,3)); 
dist_region_road_surface_region_ids = dist_region_road_surface(mask, 2);

%
region_ids = unique(dist_region_road_surface_region_ids);
color = cmap(1:10:10*numel(region_ids),:);
n = 0
for i = 1:numel(region_ids)
    
    mask = (region_ids(i) - 0.5 <= dist_region_road_surface_region_ids)&(dist_region_road_surface_region_ids <= region_ids(i) + 0.5 );
    h1 = bar(region_ids(i), sum(mask), 'barwidth', 0.75);
    if i == 1, hold on, end
    h1.FaceColor = color(i,:);
    if sum(mask) > n
        n = sum(mask)
    end

end

%
x_range = [min(region_ids)-1:1: max(region_ids) + 1];
y_interval = 50;
y_num_space = ceil(max(n)/y_interval);
y_range = [0 : y_interval: y_num_space*y_interval];

xLabel = cell(1,numel(x_range));

for j = 1:numel(x_range)
    if (j == 1 )|(j == numel(x_range))
        xLabel{j} = '';
    else
        xLabel{j} = sprintf('%.0f',x_range(j));

        
    end
end


% Set lims
ax = gca;
ax.TickLength = [0.01, 0.01];
ax.XAxis.TickLabelFormat = '%.0f';
ax.YAxis.TickLabelFormat = '%.0f';
ax.FontName = 'Time New Roman';
ax.FontSize = 24;
ax.XLim = [min(x_range), max(x_range)];
ax.YLim = [min(y_range), max(y_range)];
ax.XTick = x_range;
ax.YTick = y_range;

ax.XAxis.TickLabel = xLabel;

txt = ax.XLabel;
txt.String = 'Segment';
txt.FontSize = 24;

txt = ax.YLabel;
txt.String = 'Number of cells';
txt.FontSize = 24;

grid on

% color = cmap(1:25:64,:)
% h1 = findobj(gca,'Type','patch');
% h1.FaceColor = color;
% h1.EdgeColor = 'w';
% h1.LineWidth = 5;
file_name = sprintf('Distribution of cells within the region');
F = getframe(gcf);
imwrite(F.cdata,strcat(file_name,'.tif'),'Resolution',300)
saveas(gcf,strcat(file_name,'.png'))


%% Plot results of Ssupstr.bot
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 900 900]);
set(gca,'DataAspectRatio',[1 1 1]);
hold all
% colormap('lines');
colormap('jet');
cmap = colormap;

% Plot the region
for region_id = 5:length(SuperStructure.Component)
    comp_ptc_ids = vertcat(SuperStructure.Component(region_id).cell.ptc_ids);
    comp_ptc_xyz = OQTR.pts(comp_ptc_ids,1:3);
    
    % Plot region
    plot3(comp_ptc_xyz(:,1), comp_ptc_xyz(:,2), comp_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap((region_id-4)*12,:));
    
end

axis off
view(-70,25)
file_name = sprintf('Superstructure Bottom Surface');
F = getframe(gcf);
imwrite(F.cdata,strcat(file_name,'.tif'),'Resolution',300)
saveas(gcf,strcat(file_name,'.png'))
%% Plot intermediate surface

% Candidate points
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 900 900]);
set(gca,'DataAspectRatio',[1 1 1]);
hold all
colormap('lines');
cmap = colormap;

% Plot the points within the cell

for j = 1:numel(intermediate_surf.cell)
    cell_pts_ids = intermediate_surf.cell(j).ptc_ids;
    cell_pts_xyz = OQTR.pts(cell_pts_ids,1:3);

    color_id = j - floor(j/size(cmap,1))*size(cmap,1) + 1;
    plot3(cell_pts_xyz(:,1), cell_pts_xyz(:,2), cell_pts_xyz(:,3),'.', 'markersize', 20,'color',cmap(color_id,:));
 
end

% Plot the top and bottom
for j = 5:length(SuperStructure.Component)
    num_cells = length(SuperStructure.Component(j).cell);
%     for j = 1:num_cells
        % For each cell
        cell_ptc_ids = vertcat(SuperStructure.Component(j).cell.ptc_ids);
        cell_ptc_xyz = OQTR.pts(cell_ptc_ids,1:3);

        % Plot each cells
        color_id = j - floor(j/size(cmap,1))*size(cmap,1)+1;
        plot3(cell_ptc_xyz(:,1), cell_ptc_xyz(:,2), cell_ptc_xyz(:,3),'.', 'markersize', 1,'color', [220, 220, 220]/255);%cmap(color_id,:));
%     end
end

axis off
view(-70,25)

file_name = sprintf('Super_Bot Surf_CanPoints_Inter_Surf');
F = getframe(gcf);
imwrite(F.cdata,strcat(file_name,'.tif'),'Resolution',300)
saveas(gcf,strcat(file_name,'.png'))

%% Plot distribution of distance from cells to the road center
bw = 0.2;
[fi,zi,~] = ksdensity(dist_cells_road,'npoints',100,'bandwidth',bw,'Kernel','epanechnikov');
peak_shape_dist = peak_shape_width(zi, fi);
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 900 900]);
% set(gca,'DataAspectRatio',[1 1 1]);
hold all
% colormap('lines');
colormap('jet');
cmap = colormap;

% Plot the distance distribution

% Plot hist based distance
% hist(dist_region_road_surface(:,1), num_bin);
num_bin = ceil((max(dist_cells_road) - min(dist_cells_road))/bw);
[n, xout] = hist(dist_cells_road, num_bin);
bin_width = diff(xout);

color = cmap(ceil(linspace(1, 63, numel(n))),:);

for i = 1:numel(n)
    h1 = bar(xout(i), n(i), 'barwidth', bin_width(1) );
    if i == 1, hold on, end
    h1.FaceColor = color(i,:);

end
% h1 = bar(xout, n)
% 
% % Add KDE

plot(zi, 10*fi, '-', 'linewidth', 2,'color','b');

% Set label
% x1 = min(min(xout) - bin_width(1), min(zi));
% x2 = max(max(xout) + - bin_width(1), max(zi));
% x_num_space = ceil((x2-x1)/(0.5*bin_width(1)));
% 
% x_range = x1:bin_width(1)/2:x1+(x_num_space)*bin_width(1)/2;

xout_ = xout(n ~=0);
x_id = find(diff(xout_) >= 1);
x_id_1 = [1, x_id + 1];
x_id_2 = [x_id, numel(xout_)];
%xout_ = xout_(x_id);
xout_ = (xout_(x_id_1) + xout_(x_id_2))/2;
x_range = sort([min(zi), xout_, max(zi)]);


y_interval = 2.5;
y_num_space = ceil(max(n)/y_interval);
y_range = [0 : y_interval: y_num_space*y_interval];




% Set lims
ax = gca;
ax.TickLength = [0.01, 0.01];
ax.XAxis.TickLabelFormat = '%.2f';
ax.YAxis.TickLabelFormat = '%.0f';
ax.FontName = 'Time New Roman';
ax.FontSize = 24;
ax.XLim = [min(x_range), max(x_range)];
ax.YLim = [min(y_range), max(y_range)];
ax.XTick = sort(xout_);
ax.YTick = y_range;
ax.XGrid = 'off';
ax.YGrid = 'on';

txt = ax.XLabel;
txt.String = 'Sign distances from cells to a traffic direction (m)';
txt.FontSize = 24;

txt = ax.YLabel;
txt.String = 'Number of cells';
txt.FontSize = 24;


% color = cmap(1:25:64,:)
% h1 = findobj(gca,'Type','patch');
% h1.FaceColor = color;
% h1.EdgeColor = 'w';
% h1.LineWidth = 5;
% file_name = sprintf('Distance_bottom_surface_Road');
% F = getframe(gcf);
% imwrite(F.cdata,strcat(file_name,'.tif'),'Resolution',300)
% saveas(gcf,strcat(file_name,'.png'))




%% Plot hist based distance
count = 2;
mask = (peak_shape(count,1) - peak_shape(count,2) <= dist_region_road_surface(:,1))&...
           (dist_region_road_surface(:,1) <= peak_shape(count,1) + peak_shape(count,3)); 
dist_region_road_surface_region_ids = dist_region_road_surface(mask, 2);

%
region_ids = unique(dist_region_road_surface_region_ids);
color = cmap(1:10:10*numel(region_ids),:);
n = 0
for i = 1:numel(region_ids)
    
    mask = (region_ids(i) - 0.5 <= dist_region_road_surface_region_ids)&(dist_region_road_surface_region_ids <= region_ids(i) + 0.5 );
    h1 = bar(region_ids(i), sum(mask), 'barwidth', 0.75);
    if i == 1, hold on, end
    h1.FaceColor = color(i,:);
    if sum(mask) > n
        n = sum(mask)
    end

end

%
x_range = [min(region_ids)-1:1: max(region_ids) + 1];
y_interval = 50;
y_num_space = ceil(max(n)/y_interval);
y_range = [0 : y_interval: y_num_space*y_interval];

xLabel = cell(1,numel(x_range));

for j = 1:numel(x_range)
    if (j == 1 )|(j == numel(x_range))
        xLabel{j} = '';
    else
        xLabel{j} = sprintf('%.0f',x_range(j));

        
    end
end


% Set lims
ax = gca;
ax.TickLength = [0.01, 0.01];
ax.XAxis.TickLabelFormat = '%.0f';
ax.YAxis.TickLabelFormat = '%.0f';
ax.FontName = 'Time New Roman';
ax.FontSize = 24;
ax.XLim = [min(x_range), max(x_range)];
ax.YLim = [min(y_range), max(y_range)];
ax.XTick = x_range;
ax.YTick = y_range;

ax.XAxis.TickLabel = xLabel;

txt = ax.XLabel;
txt.String = 'Segment';
txt.FontSize = 24;

txt = ax.YLabel;
txt.String = 'Number of cells';
txt.FontSize = 24;

grid on

% color = cmap(1:25:64,:)
% h1 = findobj(gca,'Type','patch');
% h1.FaceColor = color;
% h1.EdgeColor = 'w';
% h1.LineWidth = 5;
file_name = sprintf('Distribution of cells within the region');
F = getframe(gcf);
imwrite(F.cdata,strcat(file_name,'.tif'),'Resolution',300)
saveas(gcf,strcat(file_name,'.png'))


























%% % Plot the the parapet

% for j = 1:length(parapet.cell)
% %     % Points of the cells
% %     parapet_cell_pts_ids = parapet.cell(j).ptc_ids;
% %     parapet_cell_pts_xyz = Tree.pts(parapet_cell_pts_ids,1:3);
% %     
% %     % Plot points in cells
% %     plot3(parapet_cell_pts_xyz(:,1), parapet_cell_pts_xyz(:,2), parapet_cell_pts_xyz(:,3),'.', 'markersize', 5,'color',cmap(ids(j),:));
%     
%     % Plot cell centers
%     parapet_cell_bounds_cent = parapet.cell(j).surface_features;
%     plot3(parapet_cell_bounds_cent(:,1), parapet_cell_bounds_cent(:,2), max(footpath_pts_range(:,3))+0.1,'.', 'markersize', 30,'color',cmap(ids(j),:));
%  
% end
% 

% Plot final parapet
for count = 1:length(final_parapet)
    
    parapet_ptc_ids = vertcat(final_parapet(count).cell.ptc_ids);
    parapet_ptc_xyz = Tree.pts(parapet_ptc_ids,1:3);
    plot3(parapet_ptc_xyz(:,1), parapet_ptc_xyz(:,2), parapet_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(count + 5,:));
%     
end
            
            






%%






% plot3(mean(cell_ptc_xyz(:,1))-fi/10, ones(numel(zi),1)*mean(cell_ptc_xyz(:,2)), zi,'-', 'linewidth', 2,'color','r')





% 

% view(105,18)
% Plot original data

subplot(1,2,1)
hold all


% Plot the line
% plot_style = {'b-', 'k-', 'r-', 'm-', 'c-'};
patch_ptc_ids = [];
for i = 1:size(peak_shape,1)
    for j = 2:3
        if j == 2
            z1_1 = peak_shape(i,1) - peak_shape(i,j);
        else
            z1_1 = peak_shape(i,1) + peak_shape(i,j);
        end
        % find fi
        mask = kde(:,3) == z1_1;
        kde_x1 = kde(mask,1);
        kde_y1 = kde(mask,2);
        mask = st_kde(:,3) == z1_1;
        st_kde_x1 = st_kde(mask,1);
        st_kde_y1 = st_kde(mask,2);
        plot3([kde_x1; st_kde_x1],[kde_y1;st_kde_y1],[z1_1, z1_1], '-', 'linewidth', 2,'color',cmap(i+2,:))
        
%         Plot the point on KDE
        plot3(kde_x1, kde_y1, z1_1, 'marker', 'd', 'markersize', 8, 'markerface', cmap(2,:))
        
        
    end
    
    % Plot the points of the local surface
%     z1 = peak_shape(i,1) - peak_shape(i,2);
%     z2 = peak_shape(i,1) + peak_shape(i,3);
%     mask = (z1 <= cell_ptc_xyz(:,3))&(cell_ptc_xyz(:,3) <= z2);
%     plot3(cell_ptc_xyz(mask,1), cell_ptc_xyz(mask,2), cell_ptc_xyz(mask,3),'.', 'markersize', 10,'color',cmap(i+2,:))

    patch_ptc_ids = [patch_ptc_ids;find(mask)];
end
% Plot original points
plot3(cell_ptc_xyz(:,1), cell_ptc_xyz(:,2), cell_ptc_xyz(:,3),'.', 'markersize', 5,'color',[128,128,128]/255)

% Plot unsurface points
% mask = ismember([1:size(cell_ptc_xyz,1)],patch_ptc_ids);
% plot3(cell_ptc_xyz(~mask,1), cell_ptc_xyz(~mask,2), cell_ptc_xyz(~mask,3),'.', 'markersize', 3,'color',[128,128,128]/255)



% lg = legend(h1,{sprintf('KDE')},'Location','north');
lg = legend([h1, h2],[{sprintf('KDE')}, {sprintf('Second derivative of KDE')}],'Location','northwest');
lg.Orientation = 'horizontal';
lg.FontSize = 18;
lg.FontName = 'Time New Roman';
lg.NumColumns = 1;
lg.Box = 'off';
    
    
% Set lims
ax = gca;
ax.TickLength = [0.01, 0.01];
ax.XAxis.TickLabelFormat = '%.1f';
ax.YAxis.TickLabelFormat = '%.1f';
ax.ZAxis.TickLabelFormat = '%.1f';
ax.FontName = 'Time New Roman';
ax.FontSize = 18;
ax.XLim = [min(x_range), max(x_range)];
ax.YLim = [min(y_range), max(y_range)];
ax.ZLim = [min(z_range), max(z_range)];
ax = gca;
ax.XTick = x_range;
ax.YTick = y_range;
ax.ZTick = z_range;

txt = ax.XLabel;
txt.String = 'x';
txt.FontSize = 18;

txt = ax.YLabel;
txt.String = 'y';
txt.FontSize = 18;

txt = ax.ZLabel;
txt.String = 'z';
txt.FontSize = 18;


grid on
view(-115,10)





%%
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 900 900]);
set(gca,'DataAspectRatio',[0.25 1 1]);
hold all
colormap('lines');
cmap = colormap;

count = hist(dist_parapet_road,dist_bin_edges)



axis_interval = [0.1, 0.1, 0.05];
x_range = [0, max(dist_bin_edges)];
y_range = [min(cell_sub_pts_xyz(:,2)), max(cell_sub_pts_xyz(:,2))];
z_range = [min(cell_sub_pts_xyz(:,3)), max(max(cell_sub_pts_xyz(:,3)), nd_surf_cent_2(3))];

x_range = [floor(x_range(1)*10)/10, ceil(x_range(2)*10)/10];
y_range = [floor(y_range(1)*10)/10, ceil(y_range(2)*10)/10];
z_range = [floor(z_range(1)*10)/10, ceil(z_range(2)*10)/10];

xyz_cent = [mean(x_range), mean(y_range), mean(z_range)];
num_space = ([diff(x_range), diff(y_range), diff(z_range)]./axis_interval);
% num_space = ceil([diff(x_range), diff(y_range), diff(z_range)]./axis_interval);
% num_space = 2*ceil(num_space/2);

x_range = xyz_cent(1)-axis_interval(1)*num_space(1)/2:axis_interval(1):xyz_cent(1)+axis_interval(1)*num_space(1)/2;
y_range = xyz_cent(2)-axis_interval(2)*num_space(2)/2:axis_interval(2):xyz_cent(2)+axis_interval(2)*num_space(2)/2;
z_range = xyz_cent(3) -axis_interval(3)*num_space(3)/2:axis_interval(3):xyz_cent(3)+axis_interval(3)*num_space(3)/2;




%% Plot extracted surface
subplot(1,2,2)
hold all
for j = 1:size(peak_shape,1)
    % Extract the points within the shape
    mask = (peak_shape(j,1) - peak_shape(j,2)/2 <= cell_ptc_xyz(:,3))&(cell_ptc_xyz(:,3) <= peak_shape(j,1) + peak_shape(j,2)/2);
    peak_ptc_ids = cell_ptc_ids(mask);
    peak_ptc_xyz = cell_ptc_xyz(mask,:);
    if numel(peak_ptc_ids) >= threshold.min_num_pts
        plot3(peak_ptc_xyz(:,1), peak_ptc_xyz(:,2), peak_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(j,:))
    end

end  
view(-20,15)    
%%
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 900 900]);
set(gca,'DataAspectRatio',[1 1 1]);
hold on
colormap('lines');
cmap = colormap;
plot3(peak_ptc_xyz(:,1), peak_ptc_xyz(:,2), peak_ptc_xyz(:,3),'.', 'markersize', 10,'color',cmap(j,:))
 



if numel(peak_ptc_ids) >= threshold.min_num_pts
    % Filter the surface
    [w_ptc_center, w_ptc_normal] = wPCA('ptc', peak_ptc_xyz, 'angle_threshold', 0.1*pi/180, 'ini_normal_vector', []);
    dist_ptc_plane = dist_3Dpoints_3Dplane(peak_ptc_xyz, [w_ptc_center, w_ptc_normal]);
    mask  = abs(dist_ptc_plane) <= threshold.max_distance;
    plane_ptc_ids = peak_ptc_ids(mask);
    plane_ptc_xyz = peak_ptc_xyz(mask, :);


    
end
        
plot3(plane_ptc_xyz(:,1), plane_ptc_xyz(:,2), plane_ptc_xyz(:,3),'.', 'markersize', 20,'color','r')

plot3(w_ptc_center(:,1), w_ptc_center(:,2), w_ptc_center(:,3),'.', 'markersize', 20,'color','b')



%%
start_col_ptc_ids = vertcat(Column(start_col).surface.ptc_ids);
start_col_surf_ptc_ids = Column(start_col).surface(start_col_surf).ptc_ids;
start_col_ptc_ids = setdiff(start_col_ptc_ids, start_col_surf_ptc_ids);
end_col_ptc_ids = vertcat(Column(end_col).surface.ptc_ids);
end_col_surf_ptc_ids = Column(end_col).surface(end_col_surf).ptc_ids;
end_col_ptc_ids = setdiff(end_col_ptc_ids, end_col_surf_ptc_ids);

start_col_ptc_xyz = ptc(start_col_ptc_ids,:);
start_col_surf_ptc_xyz = ptc(start_col_surf_ptc_ids,:);
end_col_ptc_xyz = ptc(end_col_ptc_ids,:);
end_col_surf_ptc_xyz = ptc(end_col_surf_ptc_ids,:);

%% Plot the columns

h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 900 900]);
set(gca,'DataAspectRatio',[1 1 1]);
hold all
colormap('lines');
cmap = colormap;

% Plot column 1
plot3(start_col_ptc_xyz(:,1), start_col_ptc_xyz(:,2), start_col_ptc_xyz(:,3),'.', 'markersize', 3,'color',[128,128,128]/255)

plot3(end_col_ptc_xyz(:,1), end_col_ptc_xyz(:,2), end_col_ptc_xyz(:,3),'.', 'markersize', 3,'color',[128,128,128]/255)

plot3(start_col_surf_ptc_xyz(:,1), start_col_surf_ptc_xyz(:,2), start_col_surf_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(1,:))

plot3(end_col_surf_ptc_xyz(:,1), end_col_surf_ptc_xyz(:,2), end_col_surf_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(3,:))

% Plot candidate beam points
% plot3(beam_filter_ptc_xyz(:,1), beam_filter_ptc_xyz(:,2), beam_filter_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(5,:))


% Plot column 1
offset = 4.0;
plot3(start_col_ptc_xyz(:,1), start_col_ptc_xyz(:,2)-offset, start_col_ptc_xyz(:,3),'.', 'markersize', 3,'color',[128,128,128]/255)

plot3(end_col_ptc_xyz(:,1), end_col_ptc_xyz(:,2)-offset, end_col_ptc_xyz(:,3),'.', 'markersize', 3,'color',[128,128,128]/255)

plot3(start_col_surf_ptc_xyz(:,1), start_col_surf_ptc_xyz(:,2)-offset, start_col_surf_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(1,:))

plot3(end_col_surf_ptc_xyz(:,1), end_col_surf_ptc_xyz(:,2)-offset, end_col_surf_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(3,:))

% Plot candidate beam points
plot3(beam_filter_ptc_xyz(:,1), beam_filter_ptc_xyz(:,2)-offset, beam_filter_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(5,:))


axis off

view(-19,50)
F = getframe(gcf);
imwrite(F.cdata,'Column_beam.tif','Resolution',300)

% Write data for future
outputFile = 'start_colum.txt';
write_txt(outputFile, start_col_ptc_xyz)
outputFile = 'start_colum_surface.txt';
write_txt(outputFile, start_col_surf_ptc_xyz)
outputFile = 'end_colum.txt';
write_txt(outputFile, end_col_ptc_xyz)
outputFile = 'end_colum_surface.txt';
write_txt(outputFile, end_col_surf_ptc_xyz)

%% Plot the beam + segmentation
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 900 900]);
set(gca,'DataAspectRatio',[1 1 1]);
hold all
colormap('lines');
cmap = colormap;
% Plot beam
% beam_segment = ptc_segment_info;
offset = 0.;
for i = 1:max(beam_segment(:,2))
    mask = beam_segment(:,2) == i;
    segment_ptc_ids = beam_segment(mask,1);
    segment_ptc_xyz = beam_filter_ptc_xyz(segment_ptc_ids,:);
    plot3(segment_ptc_xyz(:,1), segment_ptc_xyz(:,2)-offset, segment_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(i,:))
end

% Plot segment conectivity
% beam_connectivity = ptc_segment_info;
offset = 2.0;
for i = 1:max(beam_connectivity(:,2))
    mask = beam_connectivity(:,2) == i;
    segment_ptc_ids = beam_connectivity(mask,1);
    segment_ptc_xyz = beam_filter_ptc_xyz(segment_ptc_ids,:);
    plot3(segment_ptc_xyz(:,1), segment_ptc_xyz(:,2)-offset, segment_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(i,:))
end

% Plot segment conectivity
% beam_filter = ptc_segment_info;
offset = 4.0;
for i = 1:max(beam_filter(:,2))
    mask = beam_filter(:,2) == i;
    segment_ptc_ids = beam_filter(mask,1);
    segment_ptc_xyz = beam_filter_ptc_xyz(segment_ptc_ids,:);
    plot3(segment_ptc_xyz(:,1), segment_ptc_xyz(:,2)-offset, segment_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(i,:))
end

axis off

view(-19,50)
F = getframe(gcf);
imwrite(F.cdata,'Column_beam_filter.tif','Resolution',300)

% Write data for future
outputFile = 'beam_cadidate_points.txt';
write_txt(outputFile, beam_filter_ptc_xyz)
outputFile = 'beam_segments_after_voxel_region.txt';
write_txt(outputFile, [beam_filter_ptc_xyz(beam_segment(:,1),:), beam_segment(:,2)]);
outputFile = 'beam_segments_after_surface_connectivity.txt';
write_txt(outputFile, [beam_filter_ptc_xyz(beam_connectivity(:,1),:), beam_connectivity(:,2)]);
outputFile = 'beam_segments_after_filtering.txt';
write_txt(outputFile, [beam_filter_ptc_xyz(beam_filter(:,1),:), beam_filter(:,2)]);


%% Plot columns
leaf_cell_ids = Node_Leaf(Tree);
comp_cell_ids = comp_region_cell_info(mask,1);

% Neighbour
neighbour_cell_ids = Query_Neighbour_Cells(Tree,comp_cell_ids,leaf_cell_ids);
        

comp_ptc_ids = vertcat(Tree.cell_pts(comp_cell_ids).id);
comp_ptc_xyz = Tree.pts(comp_ptc_ids,1:3);
%%
h = get(0,'Children');
if isempty(h)
    figure(1)
else
    figure(max(vertcat(h.Number))+1);
end
fig = gcf;
set(fig,'Color','w');
set(fig,'Units', 'pixels');
set(fig,'Renderer','Zbuffer');
set(fig,'OuterPosition',[50 50 1200 1000]);
set(gca,'DataAspectRatio',[1 1 1]);
hold all
colormap('lines');
cmap = colormap;
view_angle = 30;
% Plot the cells around the column
% col_id = [];
% for i = 1:numel(neighbour_cell_ids)
%     if any(ismember(comp_cell_ids,neighbour_cell_ids(i)) )
%         col_id = [col_id;i];
%     end
%     neighbour_cell_ptc_ids = Tree.cell_pts(neighbour_cell_ids(i)).id;
%     neighbour_cell_ptc_xyz = Tree.pts(neighbour_cell_ptc_ids,1:3);
%     plot3(neighbour_cell_ptc_xyz(:,1), neighbour_cell_ptc_xyz(:,2), neighbour_cell_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(i,:));
% end

% % Plot the candidate points of the column
% offset_x = 3;
% offset_y = tan(deg2rad(view_angle))*offset_x;1.5;
% for i = 1:numel(comp_cell_ids)
%     comp_cell_ptc_ids = Tree.cell_pts(comp_cell_ids(i)).id;
%     comp_cell_ptc_xyz = Tree.pts(comp_cell_ptc_ids,1:3);
%     plot3(comp_cell_ptc_xyz(:,1)+offset_x, comp_cell_ptc_xyz(:,2)+offset_y, comp_cell_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(col_id(i),:));
% end

% % Plot segment
% % column_segment = ptc_segment_info;
% offset_x = 6;
% offset_y = tan(deg2rad(view_angle))*offset_x;1.5;
% for i = 1:max(column_segment(:,2))
%     mask = column_segment(:,2) == i;
%     segment_ptc_ids = column_segment(mask,1);
%     segment_ptc_xyz = comp_ptc_xyz(segment_ptc_ids,:);
%     plot3(segment_ptc_xyz(:,1)+offset_x, segment_ptc_xyz(:,2)+offset_y, segment_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(i,:))
% end
% 
% column_connectivity = ptc_segment_info;
% offset_x = 9;
% offset_y = tan(deg2rad(view_angle))*offset_x;1.5;
% for i = 1:max(column_connectivity(:,2))
%     if any(i==[1,4])
%     mask = column_connectivity(:,2) == i;
%     segment_ptc_ids = column_connectivity(mask,1);
%     segment_ptc_xyz = comp_ptc_xyz(segment_ptc_ids,:);
%     plot3(segment_ptc_xyz(:,1)+offset_x, segment_ptc_xyz(:,2)+offset_y, segment_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(i,:))
%     end
% end
% 
% % column_outlier_filter = ptc_segment_info;
offset_x = 12;
offset_y = tan(deg2rad(view_angle))*offset_x;1.5;
for i = 1:max(column_outlier_filter(:,2))
    if any(i==[1,4])
    mask = column_outlier_filter(:,2) == i;
    segment_ptc_ids = column_outlier_filter(mask,1);
    segment_ptc_xyz = comp_ptc_xyz(segment_ptc_ids,:);
    plot3(segment_ptc_xyz(:,1)+offset_x, segment_ptc_xyz(:,2)+offset_y, segment_ptc_xyz(:,3),'.', 'markersize', 5,'color',cmap(i,:))
    end
end
axis off
view(view_angle,60)
% %
% F = getframe(gcf);
% imwrite(F.cdata,'Column_Etraction_Surface_Filter.tif','Resolution',300)

%%
% plot histograme
format bank
numBinX = 4;
BinWidth = (max(sign_dist)-min(sign_dist))/numBinX;
h = gcf;
if h == 1
    figure(h)
else
    figure(h+1)
end
histogram(sign_dist,numBinX)
[n,xout] = hist(sign_dist,numBinX);
ax = gca;
ax.FontSize = 16;
% set(gca,'YLim',[0 40]);

% set(gca,'XLim',[min(xout)-0.1*BinWidth max(xout)+0.1*BinWidth]);
set(gca,'XLim',[min(xout)-1*BinWidth, max(xout)+1*BinWidth]);
set(gca,'XTick',min(xout)-1*BinWidth:BinWidth:max(xout)+1*BinWidth);
% NumXTick = ((max(xout))-(min(xout)))/(0.5*BinWidth);
%
xTickMarks = cell(1,NumXTick);
for i=1:(NumXTick+2)
    if i == 2
        xTickMarks(i) = {-1};
    elseif i == 5
        xTickMarks(i) = {1};
    else
        xTickMarks(i) = {''};
    end
%     if any(i ==[2,5])
% %         temp = round((min(xout)-1.0*BinWidth+(i-1)*BinWidth)*10)
% %         xTickMarks(i) = {num2str(temp/10,2)};
%         xTickMarks(i) = {round(
% %     else
% %         xTickMarks(i) = {''}; 
%     end
end
set(gca,'XTickLabel',xTickMarks);
txt = ax.XLabel;
txt.String = 'Sign of Distance';
    


