function peak_shape = peak_shape_width(x, y)
%% This function is to estimate the width of the shape of the local peaks. 
% This as done by using Gausian to fit the points within the peak shape
% This function was used some functions from http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm 
% The end of the peak where the slope change significantly, which means the
% second derivative is maximum
% Input:
% Output:
% Developed:
% 
% Demo:

    % x = zi;
    % y = fi;

%%
% Determine a number of peaks from KDE
[~, pks_locs] = findpeaks(y);

% Compute the second derivative
y_deriv = deriv(x, y);
y_deriv_2 = deriv(x, y_deriv);

% Preallocation
peak_shape = zeros(numel(pks_locs),3);
peak_shape(:,1) = x(pks_locs);

for i = 1:numel(pks_locs)
    % Determine a location of the peak
    pk_loc = pks_locs(i);
    
    % Searching points on the left
    left_ids = 1:pk_loc-1;
    mask = (y_deriv(left_ids) <= 0);
    if any(mask)
        start_peak = left_ids(find(mask==1, 1, 'last'));
    else
        start_peak = left_ids(1);
    end
    x_left = x(start_peak:pk_loc-1);
%     y_left = y(start_peak:pk_loc);
    y_left_deriv2 = y_deriv_2(start_peak:pk_loc-1);
    [~, mask] = max(y_left_deriv2);
    x_left = x_left(mask);
    
    % Searching points on the right
    right_ids = pk_loc+1:numel(y);
    mask = (y_deriv(right_ids) >= 0);
    if any(mask)
        end_peak = right_ids(find(mask==1, 1, 'first'));
    else
        end_peak = right_ids(end);
    end
    x_right = x(pk_loc+1:end_peak);
%     y_right = y(pk_loc:end_peak);
    y_right_deriv2 = y_deriv_2(pk_loc+1:end_peak);
    [~, mask] = max(y_right_deriv2);
    x_right = x_right(mask);
    

    % Estimate width of the peaks by gausian fitting
%     if (numel(x_peak) > 3) && (numel(y_peak) > 3)
% %         [~, ~, peak_width] = gaussfit(x_peak,y_peak);
%         peak_width = max(x_peak) - min(x_peak);
%     else
%         peak_width = 0;
%     end
%     peak_width = abs(x(start_peak) - x(end_peak));
%     peak_shape(i,2) = x_right - x_left;
    peak_shape(i,[2,3]) = [x(pk_loc) - x_left, x_right - x(pk_loc)];

end

% 
% function [Height, Position, Width] = gaussfit(x,y)
% % Converts y-axis to a log scale, fits a parabola
% % (quadratic) to the (x,ln(y)) data, then calculates
% % the position, width, and height of the
% % Gaussian from the three coefficients of the
% % quadratic fit.  This is accurate only if the data have
% % no baseline offset (that is, trends to zero far off the
% % peak) and if there are no zeros or negative values in y.
% %
% % Example 1: Simplest Gaussian data set
% % [Height, Position, Width]=gaussfit([1 2 3],[1 2 1]) 
% %    returns Height = 2, Position = 2, Width = 2
% %
% % Example 2: best fit to synthetic noisy Gaussian
% % x=50:150;y=100.*gaussian(x,100,100)+10.*randn(size(x));
% % [Height,Position,Width]=gaussfit(x,y) 
% %   returns [Height,Position,Width] clustered around 100,100,100.
% %
% % Example 3: plots data set as points and best-fit Gaussian as line
% % x=[1 2 3 4 5];y=[1 2 2.5 2 1];
% % [Height,Position,Width]=gaussfit(x,y);
% % plot(x,y,'o',linspace(0,8),Height.*gaussian(linspace(0,8),Position,Width))
% % Copyright (c) 2012, Thomas C. O'Haver
% 
%     maxy = max(y);
%     for p=1:length(y)
%         if y(p)<(maxy/100)
%             y(p)=maxy/100;
%         end
%     end % for p=1:length(y),
%     log_y = log(y);
%     coef = polyfit(x, log_y, 2);
%     a = coef(3);
%     b = coef(2);
%     c = coef(1);
%     Height = exp(a-c*(b/(2*c))^2);
%     Position = -b/(2*c);
%     Width = 2.35482/(sqrt(2)*sqrt(-c));

% function d = deriv(x, y)
%     % First derivative of vector using 2-point central difference.
%     
%     n=length(x);
%     d(1) = (y(2) - y(1))/(x(2) - x(1));
%     d(n) = (y(n) - y(n-1))/(x(n) - x(n-1));
%     for j = 2:n-1
%         d(j)=(y(j+1)-y(j))./(x(j+1) - x(j));
% %         d(j)=(y(j+1)-y(j-1))./(x(j+1) - x(j-1));
%     end
    
% function [backward, forward] = slope(x, y)
%     % Slope at the points: backward and forward.
%     
%     n = length(x);
%     backward(1) = (y(1) - y(2))/(x(1) - x(2));
%     backward(n) = (y(n-1) - y(n))/(x(n-1) - x(n));
%     forward(1) = (y(2) - y(1))/(x(2) - x(1));
%     forward(n) = (y(n) - y(n-1))/(x(n) - x(n-1));
%     for j = 2:n-1
%         backward(j)=(y(j-1)-y(j))./ (x(j-1) - x(j));
%         forward(j)=(y(j+1)-y(j))./ (x(j+1) - x(j));
%     end

% function ratio = slope_ratio(y_slope)
%     % Compute a ratio of the slope at the points 
%     n = length(y_slope);
%     ratio(1) = 1;
%     ratio(n) = 1;
%     ratio(2:n-1) = y_slope(1:end-2,:)./y_slope(3:end,:);%max([y_slope(3:end,:)./y_slope(1:end-2,:),y_slope(1:end-2,:)./y_slope(3:end,:)],[],2);
%     
%%
    
%     figure(2)
%     hold all
%     plot(x,y,'r-')
%     plot(x,y_deriv,'b-')
%     plot(x,y_deriv_2,'m-')
    
        
