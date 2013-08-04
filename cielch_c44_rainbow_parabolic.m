% Colour difference is not the same throughout
% Colours are not saturated enough (chroma too low)

% First can be readily fixed, second cannot

function rgb = cielch_c44_rainbow_parabolic(n, func, debug)

% Input handling
if nargin<1
    n = size(get(gcf,'colormap'),1);
end
if nargin<2
    func = [];
end
if nargin<3
    debug = 0;
end

use_uplab = false;

chr = 44; % Lhc chroma of 44

% Curve co-ordinates in Lch space for constant chroma
xx = [0 152.5 305];
yy = [22 90 22];
% xx = [35, 167.5, 300];
% yy = [16.5, 90, 16.5];
% xx = [0 152.5 305];
% yy = [40 82 40];


% Parabola center
x0 = xx(2); % 167.5;
y0 = yy(2); % 91;

% (y-y0) = p4*(x-x0).^2

p4 = (yy(1)-y0) / (xx(1)-x0).^2;
f  = abs(1/(4*p4)); % focal length

h = (x0-xx(1))/2;
q = sqrt(f.^2+h.^2);
s = h*q/f + f * log((h+q)/f); % arc length of half the parabola

% s = x/2*sqrt(f.^2+x.^2/4) + f*log((x/2+sqrt(f.^2+x.^2/4))/f)
% exp(s) = exp(x/2*sqrt(f.^2+x.^2/4)) * ((x/2+sqrt(f.^2+x.^2/4))/f).^f

% dybydx = 2*p4*(x-x0)
% ds = sqrt(1+dybydx.^2)*dx;

% S = solve(...
%     '170 = x/2*sqrt(f^2+x^2/4) + f*log((x/2+sqrt(f^2+x^2/4))/f)',...
%     sprintf('f=%f',f));
% 
% S = solve(...
%     sprintf('f=%f',f),...
%     'q = sqrt(f^2+h^2)',...
%     '170 = h*q/f + f * log((h+q)/f)'...
%     );



% Want to make the points equally spaced from one another in L*-h_{a*b*}
% space. Do this by picking points based on their arc length from the
% centre of the parabola.

neach = ceil(n/2);

switch mod(n,2)
    case 0
        % n is even
        sss = linspace(0,s,n);
        sss = sss(1:2:end);
%         sss = linspace(-s,s,n);
%         sss = sss(neach+1:end)
    case 1
        % n is odd
        sss = linspace(0,s,neach);
end

% Only solve for half the parabola, as we can use the symmetry to find both
% at once
hhh = nan(size(sss));
for i=1:neach
    s_ = sss(i);
    hhh(i) = solve(...
        sprintf('%2$f = h*sqrt(%1$f^2+h^2)/%1$f + %1$f * log((h+sqrt(%1$f^2+h^2))/%1$f)',...
            f,s_)...
        );
end

xxx1 = x0 - hhh*2;
xxx1 = max(xxx1,0);
xxx1 = fliplr(xxx1);

xxx2 = hhh*2 + x0;

switch mod(n,2)
    case 0
        % n is even
%         xxx1 = xxx1(1:end-1);
%         xxx2 = xxx2(2:end  );
    case 1
        % n is odd
        xxx1 = xxx1(1:end-1);
%         xxx2 = xxx2(2:end  );
end

xxx = [xxx1 xxx2];
% Swap so that blue-red not red-blue
xxx = xxx(end:-1:1);

yyy = y0 + p4*(xxx-x0).^2;

% Move into Lhc
L = yyy;
h = xxx;
c = chr;

% Move from Lhc into Lab
a = c.*cos(h/360*(2*pi));
b = c.*sin(h/360*(2*pi));

Lab = [L' a' b'];

% Move from Lab into rgb
rgb = gd_lab2rgb(Lab, use_uplab, func);

% -- Plot for debug ---
if debug
v  = round(chr);

rgbgamut = fetch_cielchab_gamut();

I = rgbgamut.lch_chr.Lmin.lch(:,2)==v;

figure;
set(gca,'Color',[0.4663 0.4663 0.4663]);
set(gca,'XLim',[0 360]);
set(gca,'YLim',[0 100]);
xlabel('Hue (Lch_{ab})');
ylabel('Lightness (Lch_{ab})');
title(sprintf('chroma %.0f',v));
set(gca,'NextPlot','add');

% plot(...
%     [rgbgamut.lch(I1,2);-flipud(rgbgamut.lch(I2,2))],...
%     [rgbgamut.lch(I1,1); flipud(rgbgamut.lch(I2,1))],...
%     'k');
% set(main_axes,'NextPlot','add');

scatter(...
    rgbgamut.lch_chr.Lmin.lch(I,3),...
    rgbgamut.lch_chr.Lmin.lch(I,1),...
    20,...
    rgbgamut.lch_chr.Lmin.rgb(I,:),...
    'fill');
scatter(...
    rgbgamut.lch_chr.Lmax.lch(I,3),...
    rgbgamut.lch_chr.Lmax.lch(I,1),...
    20,...
    rgbgamut.lch_chr.Lmax.rgb(I,:),...
    'fill');
% plot( 35  , 16.5,'xk')
% plot(110  , 75.9,'xk')
% plot(167.5, 91  ,'xk')
% plot(225  , 75.9,'xk')
% plot(300  , 16.5,'xk')

% plot( [35, 110, 167.5, 225, 300], [16.5, 75.9, 91, 75.4, 16.5], 'ko-');
% plot( [35, 110, 225, 300], [16.5, 90, 75.4, 16.5], 'ko-');

% % plot(xxx,yyy,'-k');
scatter(...
    xxx,...
    yyy,...
    20,...
    rgb,...
    'fill');

img = repmat(rgb,[1 1 20]);
img = permute(img,[1 3 2]);

figure;
imagesc(img);
axis xy;
% figure;
% imagesc(img(1:4:end,:,:));
% axis xy;

end

end