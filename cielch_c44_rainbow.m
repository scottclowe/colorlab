function cielch_c44_rainbow(n, func)

if nargin<1
    n = size(get(gcf,'colormap'),1);
end
if nargin<2
    func = [];
end

use_uplab = false;


xx = [35, 110, 225, 300]
yy = [16.5, 90, 75.4, 16.5]
mm = [1.0375, 0, -2.6, -1.0375]

% Make curve in 3 parts, fitting four parameters each time
a = nan(4,3);
% Fit a thrid order quadratic
for i=1:3
%    y(x) = a(1,i) + a(2,i)*x + a(3,i)*x^2 + a(4,i)*x^3;
%    m(x) = a(2,i) + 2*a(3,i)*x + 3*a(4,i)*x^2;
%    
%    mm(j)         = a(2,i) + 2*a(3,i)*xx(j) + 3*a(4,i)*xx(j)^2;
%    mm(j)-mm(j+1) = 2*a(3,i)*(xx(j)-xx(j+1)) + 3*a(4,i)*(xx(j)^2-xx(j+1)^2);
%    
%    yy(j)         = a(1,i) + a(2,i)*xx(j) + a(3,i)*xx(j)^2 + a(4,i)*xx(j)^3;
%    yy(j)-yy(j+1) = a(2,i)*xx(j) + a(3,i)*xx(j)^2 + a(4,i)*xx(j)^3;
%    
%    mm(j)*xx(j)   = a(2,i)*xx(j) + 2*a(3,i)*xx(j)^2 + 3*a(4,i)*xx(j)^3;
%    3*yy(j)  -mm(j)  *xx(j)   = 3*a(1,i) + 2*a(2,i)*xx(j)   + a(3,i)*xx(j)^2;
%    3*yy(j+1)-mm(j+1)*xx(j+1) = 3*a(1,i) + 2*a(2,i)*xx(j+1) + a(3,i)*xx(j+1)^2;
%    xx(j+1)*(3*yy(j)  -mm(j)  *xx(j)) - xx(j)*(3*yy(j+1)-mm(j+1)*xx(j+1)) = a(3,i) (xx(j)^2*xx(j+1) - xx(j+1)^2*xx(j));
%    
%    a(3,i) = (xx(j+1)*(3*yy(j)  -mm(j)  *xx(j)) - xx(j)*(3*yy(j+1)-mm(j+1)*xx(j+1))) / (xx(j)^2*xx(j+1) - xx(j+1)^2*xx(j));
%    a(2,i) = (3*yy(j)  -mm(j)  *xx(j) - a(3,i)*xx(j)^2) / (2*xx(j));
%    a(4,i) = (mm(j) - (a(2,i)+ 2*a(3,i)*xx(j))) / (3*xx(j)^2);

solve(...
    sprintf('%f = a + %f*b + %f*c + %f*d', yy(1), xx(1), xx(1)^2, xx(1)^3) ,...
    sprintf('%f = a + %f*b + %f*c + %f*d', yy(2), xx(2), xx(2)^2, xx(2)^3) ,...
    sprintf('%f = b + %f*c + %f*d',        mm(1), 2*xx(1), 3*xx(1)^2) ,...
    sprintf('%f = b + %f*c + %f*d',        mm(2), 2*xx(2), 3*xx(2)^2) ...
    );
end

cmap = gd_lab2rgb(Lab, use_uplab, func);

end