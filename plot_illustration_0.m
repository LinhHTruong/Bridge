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
    


