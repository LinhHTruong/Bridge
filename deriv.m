function d = deriv(x, y)
% First derivative of vector using 2-point central difference.

n=length(x);
d(1) = (y(2) - y(1))/(x(2) - x(1));
d(n) = (y(n) - y(n-1))/(x(n) - x(n-1));
for j = 2:n-1
    d(j)=(y(j+1)-y(j))./(x(j+1) - x(j));
%         d(j)=(y(j+1)-y(j-1))./(x(j+1) - x(j-1));
end
end