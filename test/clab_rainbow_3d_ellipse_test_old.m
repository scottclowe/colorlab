
close all;

use_uplab = false;


%%
g = fetch_cielchab_gamut();

hs = unique(g.lch(:,3));
g.lch_chproj.lch = nan(length(hs),3);
g.lch_chproj.lab = nan(length(hs),3);
g.lch_chproj.rgb = nan(length(hs),3);
for i=1:length(hs)
    h = hs(i);
    li = g.lch(:,3)==h;
    indx = find(li);
    [c,I] = max(g.lch(li,2));
    g.lch_chproj.lch(i,:) = g.lch(indx(I),:);
    g.lch_chproj.lab(i,:) = g.lab(indx(I),:);
    g.lch_chproj.rgb(i,:) = g.rgb(indx(I),:);
end

%%

%          [ L     c   h ]   [L         a         b     ]
Red      = [54.75 106  41]; % 54.2908   80.8211   69.9023
Orange   = [70     86  62];
Yellow   = [97.5   94  99];
LgtGreen = [88    134 113];
MidGreen = [57.25  80 134];
Cyan     = [70     40 220]; % approx
LgtBlue  = [56.5   65 270];
DrkBlue  = [29.5  131 301];
Purple   = [56.75  88 349];
DrkPurple= [40.25  83 327];


P = [];
P(size(P,1)+1,:) = Red;
P(size(P,1)+1,:) = Orange;
P(size(P,1)+1,:) = Yellow;
P(size(P,1)+1,:) = LgtGreen;
P(size(P,1)+1,:) = MidGreen;
P(size(P,1)+1,:) = Cyan;
P(size(P,1)+1,:) = LgtBlue;
P(size(P,1)+1,:) = DrkBlue;
P(size(P,1)+1,:) = Purple;
P(size(P,1)+1,:) = DrkPurple;
P(size(P,1)+1,:) = P(1,:);


a = P(:,2).*cos(P(:,3)/360*(2*pi));
b = P(:,2).*sin(P(:,3)/360*(2*pi));

P2 = [P(:,1) a b];

rgb = soft_lab2rgb(P2, use_uplab);

figure; hold on;
set(gca,'Color',[0.4663 0.4663 0.4663]);
plot3(P2(:,2),P2(:,3),P2(:,1),'-k')
scatter3(P2(:,2),P2(:,3),P2(:,1),50,rgb);

%%

% close;
figure; hold on;
set(gca,'Color',[0.4663 0.4663 0.4663]);
scatter(...
    g.lch_chproj.lab(:,2),...
    g.lch_chproj.lab(:,3),...
    20,...
    g.lch_chproj.rgb,...
    'filled');

% cx  = 10;
% cy  = 10;
% axa = 90;
% axb = 60;
% theta = 50;
% 
% cx  = 15;
% cy  = 10;
% axa = 80;
% axb = 60;
% theta = 50;
% 
% cx  = 20;
% cy  = 0;
% axa = 95;
% axb = 55;
% theta = 60;

% C = [15 5];
% axa = 90;
% axb = 60;
% theta = 55;

C = [20 10 65];
a = 80;
b = 55;
theta = -55/180*pi;

[X Y] = calculateEllipse(C(1), C(2), a, b, theta/pi*180);
plot(X,Y,'k');
axis equal

%%
c = (X.^2 + Y.^2).^(1/2);
h = atan(Y./X); % In radians
h(isnan(h)) = 0;
h = h./(2*pi)*360; % In degrees
s = sign(X);
s(s==1) = 0;
h = h-180*s;
h = mod(h,360);

chrs = c;
max_c = max(c);
hues = h;

Lmins = nan(1,length(hues));
Lmaxs = nan(1,length(hues));
for i=1:length(hues)
    ih = round(hues(i));
    ic = round(chrs(i));
    Lmin = g.lch_chr.Lmin.lch(g.lch_chr.Lmin.lch(:,3)==ih & g.lch_chr.Lmin.lch(:,2)==ic, 1);
    Lmax = g.lch_chr.Lmax.lch(g.lch_chr.Lmax.lch(:,3)==ih & g.lch_chr.Lmin.lch(:,2)==ic, 1);
    if isempty(Lmin)
        warning('No valid L for h=%.3f c=%.3f',ih,ic);
        continue;
    end
    Lmins(i) = Lmin;
    Lmaxs(i) = Lmax;
end

hues2 = mod(hues+60,360)-60;
figure; hold on;
plot(hues2,Lmins,'k');
plot(hues2,Lmaxs,'r');
