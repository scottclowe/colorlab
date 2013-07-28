
close all;

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

cform = makecform('lab2srgb');
rgb = applycform(P2, cform);

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

%%

close all;

% Set#1: full 3d ellipse
% C = [15 10 65];
% a = 70;
% b = 50;
% theta = -70/180*pi; %atan(-(93+113)/(15+55)); % -60/180*pi;
% % phi = 0; %atan(47/180);
% phi = -atan(50/2/a);
% psi = -atan(30/2/b); %pi/4;
% U = [ cos(theta) sin(theta) sin(phi)];
% V = [-sin(theta) cos(theta) sin(psi)];

% Nearly there with partial ellipse
% C = [140 -50 15];
% a = 200;
% b = 90;
% theta = -27/180*pi; %atan(-(93+113)/(15+55)); % -60/180*pi;
% % phi = 0; %atan(47/180);
% phi = -atan((90-30)/sqrt((51+64)^2+(110-12)^2)); % -20/180*pi; % -atan(50/2/a);
% psi = 0; %5/180*pi; % atan(30/2/b); %pi/4;
% U = [ cos(theta) sin(theta) sin(phi)];
% V = [-sin(theta) cos(theta) sin(psi)];

% Working partial ellipse
% a = 300;
% b = 85;
% theta = -31/180*pi; %atan(-(93+113)/(15+55)); % -60/180*pi;
% % phi = 0; %atan(47/180);
% phi = -atan((90-30)/sqrt((51+64)^2+(110-12)^2)); % -20/180*pi; % -atan(50/2/a);
% psi = -2/180*pi; %5/180*pi; % atan(30/2/b); %pi/4;
% U = [ cos(theta) sin(theta) sin(phi)];
% V = [-sin(theta) cos(theta) sin(psi)];
% 
% % Define centre such that when t=pi, we are at P
% P = [-48 51 92];
% t1 = pi;
% C = [...
%     P(1) - (a * cos(t1) * U(1) + b * sin(t1) * V(1)) ...
%     P(2) - (a * cos(t1) * U(2) + b * sin(t1) * V(2)) ...
%     P(3) - (a * cos(t1) * U(3) + b * sin(t1) * V(3)) ...
%     ];
% t_start = 3.75;
% t_end = 1.75;

% Pretty much there ellipse
% a = 300;
% b = 94;
% theta = -30/180*pi; %atan(-(93+113)/(15+55)); % -60/180*pi;
% % phi = 0; %atan(47/180);
% phi = -atan((90-30)/sqrt((51+64)^2+(110-12)^2)); % -20/180*pi; % -atan(50/2/a);
% psi = -2/180*pi; %5/180*pi; % atan(30/2/b); %pi/4;
% t_start = 3.75;
% t_end = 1.75;
% P = [-48 57 92];

% Working ellipse
% a = 300;
% b = 91;
% theta = -25/180*pi; %atan(-(93+113)/(15+55)); % -60/180*pi;
% % phi = 0; %atan(47/180);
% phi = -17.5/180*pi; % -atan((90-30)/sqrt((51+64)^2+(110-12)^2)); % -20/180*pi; % -atan(50/2/a);
% psi = 1.7/180*pi; %5/180*pi; % atan(30/2/b); %pi/4;
% t_start = 4.31;
% t_end = 2.20;
% P = [-60 49 88];
% 
% Working under less accurate gamut
% a = 300;
% b = 93.5;
% theta = -25/180*pi; %atan(-(93+113)/(15+55)); % -60/180*pi;
% % phi = 0; %atan(47/180);
% phi = -17.5/180*pi; % -atan((90-30)/sqrt((51+64)^2+(110-12)^2)); % -20/180*pi; % -atan(50/2/a);
% psi = 1.7/180*pi; %5/180*pi; % atan(30/2/b); %pi/4;
% t_start = 4.31;
% t_end = 2.20;
% P = [-60 49 88];
% %   [  x  y  z ]
% %   [  a  b  L ]

% Working under more accurate gamut & included in cmap builder
a = 400;
b = 106;
theta = -25.5/180*pi; %atan(-(93+113)/(15+55)); % -60/180*pi;
% phi = 0; %atan(47/180);
phi = -17.5/180*pi; % -atan((90-30)/sqrt((51+64)^2+(110-12)^2)); % -20/180*pi; % -atan(50/2/a);
psi = 1.7/180*pi; %5/180*pi; % atan(30/2/b); %pi/4;
t_start = 4.2;
t_end = 2.30;
P = [-60 49 88];
%   [  x  y  z ]
%   [  a  b  L ]

% U and V are the major and minor axes respectively
% theta is the angle of rotation in the x-y plane
% phi is the angle of vertical incline of the major axis
% psi is the angle of vertical incline of the minor axis
U = [ cos(theta) sin(theta) sin(phi)];
V = [-sin(theta) cos(theta) sin(psi)];

% Define centre such that when t=pi, we are at P
t1 = pi;
C = [...
    P(1) - (a * cos(t1) * U(1) + b * sin(t1) * V(1)) ...
    P(2) - (a * cos(t1) * U(2) + b * sin(t1) * V(2)) ...
    P(3) - (a * cos(t1) * U(3) + b * sin(t1) * V(3)) ...
    ];

t0 = 60/180*pi;
n = 64;
t = t0 + linspace(0,2*pi,n);
x = C(1) + a * cos(t) * U(1) + b * sin(t) * V(1);
y = C(2) + a * cos(t) * U(2) + b * sin(t) * V(2);
z = C(3) + a * cos(t) * U(3) + b * sin(t) * V(3);

% figure;
% plot3(x,y,z);
% xlabel('x');
% ylabel('y');
% zlabel('z');

f1 = figure; hold on;
set(gca,'Color',[0.4663 0.4663 0.4663]);
scatter(...
    g.lch_chproj.lab(:,2),...
    g.lch_chproj.lab(:,3),...
    20,...
    g.lch_chproj.rgb,...
    'filled');
[X Y] = calculateEllipse(C(1), C(2), a, b, theta/pi*180);
plot(X,Y,'k');
axis equal

res = 17;

f2 = figure('Position',[150 150 550 550]);
hold on;
plot3(x,y,z,'k');
scatter3(...
    g.lab([1:res:end-1 end],2), ...
    g.lab([1:res:end-1 end],3), ...
    g.lab([1:res:end-1 end],1), ...
    20, ...
    g.rgb([1:res:end-1 end],:), ...
    'filled');
set(gca,'Color',[0.4663 0.4663 0.4663]);
set(gca,'XLim',[-150 150],'YLim',[-150 150],'ZLim',[0 100]);
xlabel('a*')
ylabel('b*')
zlabel('L*')


% start_coord = [70 -109 29];
% end_coord   = [77   70 54];
% t_start = atan2( (start_coord(2)-C(2)),(start_coord(1)-C(1)) );
% t_end   = atan2( (  end_coord(2)-C(2)),(  end_coord(1)-C(1)) );

% t_start = mod(t_start,2*pi);
% t_end = mod(t_end,2*pi);


n = 64;
t = linspace(t_start,t_end,n);
x = C(1) + a * cos(t) * U(1) + b * sin(t) * V(1);
y = C(2) + a * cos(t) * U(2) + b * sin(t) * V(2);
z = C(3) + a * cos(t) * U(3) + b * sin(t) * V(3);

figure(f2);
plot3(x(1)  , y(1)  , z(1)  , 'ko');
plot3(x(end), y(end), z(end), 'ks');

figure(f1);
plot(x(1)  , y(1)  , 'ko');
plot(x(end), y(end), 'ks');

Lab = [z' x' y'];
cform = makecform('lab2srgb');
rgb = applycform(Lab, cform);

img = repmat(rgb,[1 1 20]);
img = permute(img,[1 3 2]);

figure;
imagesc(img);
axis xy;

%
% Project onto the plane of the ellipse
% Plot chr value for each of the L and h values
% This will allow me to see how to get it to fit more easily (hopefully)

[TF,P2] = isingamut(Lab,g,'Lab');

% disp(sum(TF)/length(TF));

figure;
hold on;
plot3(Lab(:,2), Lab(:,3), Lab(:,1), 'b')
plot3(P2(:,2) , P2(:,3) , P2(:,1) , 'r')
xlabel('a');
ylabel('b');
zlabel('L');
view(0,90);


%%
% Ideas for parametrising ellipse by arc length:
%
% 1 use matlab's ODE solver: ode45
%
% 2 http://www.mathworks.co.uk/matlabcentral/fileexchange/26819-ellipse-arc-length
% arcfun = @(theta,a,b) sqrt((a.*sin(theta)).^2 + (b.*cos(theta)).^2);
% quadgk(@(t) arcfun(t,Req,Rpo),-pi/2,pi/2)
%
% 3 http://arxiv.org/pdf/math/0506384.pdf
% http://en.wikipedia.org/wiki/Ellipse#Circumference
% Ramanujan's approximation
%
% NB: built in function for ellipses only gives complete integral
% http://www.mathworks.de/de/help/matlab/ref/ellipke.html
%
% 4. Generate a ton of points and use these to find the length and points
% with appropriate distances in between
% Make a high-res map once, and interpolate from this in the future

%%

% First generate a ton of points to find the end points of theta
n = 100001;
t = linspace(t_start,t_end,n);
x = C(1) + a * cos(t) * U(1) + b * sin(t) * V(1);
y = C(2) + a * cos(t) * U(2) + b * sin(t) * V(2);
z = C(3) + a * cos(t) * U(3) + b * sin(t) * V(3);

Lab = [z' x' y'];
TF = isingamut(Lab,g,'Lab');

btsp_t_srt = t(find(TF,1,'first'));
btsp_t_end = t(find(TF,1,'last'));

% Now generate a ton of points in this range
n = 100001;
t = linspace(btsp_t_srt,btsp_t_end,n);
x = C(1) + a * cos(t) * U(1) + b * sin(t) * V(1);
y = C(2) + a * cos(t) * U(2) + b * sin(t) * V(2);
z = C(3) + a * cos(t) * U(3) + b * sin(t) * V(3);

Lab = [z' x' y'];
TF = isingamut(Lab,g,'Lab');

% Check all points are now inside
fprintf('This one needs to equal %d = %d',sum(TF),length(TF));
disp(sum(TF)/length(TF))

% Find distance between each set of consequtive points
Lab_dif = diff(Lab,1,1);

% Using simple 1931 color difference metric
Lab_inv = sqrt(sum(Lab_dif.^2,2));

% % Using more up to date color difference metric
% % https://en.wikipedia.org/wiki/Color_difference#CIEDE2000
% LCh = applycform(Lab,makecform('lab2lch'));
% Lab_inv = ciede(LCh(1:end-1,:),LCh(2:end,:));
% % However, the results produced by this do not look so good
% % Could be a mistake in the calculation

% Hence find length of complete arc
arclen_full = sum(Lab_inv);

% Find 256 points equally spaced out of the 100000 we generated
n_target = 256;
arclen_inv_target = arclen_full./(n_target-1);

% Match target to the arclength NOT the euclidean length between points
% This is easier to code
Lab_inv_runtot = [0; cumsum(Lab_inv)];
arclen_targets = arclen_inv_target * (0:(n_target-1));

t_new = interp1(Lab_inv_runtot,t,arclen_targets);

t_new(1)   = btsp_t_srt;
t_new(end) = btsp_t_end;

x = C(1) + a * cos(t_new) * U(1) + b * sin(t_new) * V(1);
y = C(2) + a * cos(t_new) * U(2) + b * sin(t_new) * V(2);
z = C(3) + a * cos(t_new) * U(3) + b * sin(t_new) * V(3);

Lab = [z' x' y'];
cform = makecform('lab2srgb');
rgb = applycform(Lab, cform);

ciebow_cmap = rgb;

img = repmat(rgb,[1 1 20]);
img = permute(img,[1 3 2]);

figure;
imagesc(img);
axis xy;
figure;
imagesc(img(1:4:end,:,:));
axis xy;

[TF,P2] = isingamut(Lab,g,'Lab');

% disp(sum(TF)/length(TF))

figure;
hold on;
plot3(Lab(:,2), Lab(:,3), Lab(:,1), 'bx')
plot3(P2(:,2) , P2(:,3) , P2(:,1) , 'r-')
xlabel('a');
ylabel('b');
zlabel('L');
view(0,90);


cmap_dif = diff(Lab,1,1);
cmap_inv = sqrt(sum(cmap_dif.^2,2));

%%
% Should add an option to have same amount of red as blue
% This is a very simple method used where the curve is split in two and
% there are different Delta E values between colours in each.
% It would be better to have the rate change sigmoidally instead

% Mid point is at t=pi

% Now generate a ton of points in this range
n = 100001;
t1 = linspace(btsp_t_srt, pi, n);
x1 = C(1) + a * cos(t1) * U(1) + b * sin(t1) * V(1);
y1 = C(2) + a * cos(t1) * U(2) + b * sin(t1) * V(2);
z1 = C(3) + a * cos(t1) * U(3) + b * sin(t1) * V(3);
t2 = linspace(pi, btsp_t_end, n);
x2 = C(1) + a * cos(t2) * U(1) + b * sin(t2) * V(1);
y2 = C(2) + a * cos(t2) * U(2) + b * sin(t2) * V(2);
z2 = C(3) + a * cos(t2) * U(3) + b * sin(t2) * V(3);

Lab1 = [z1' x1' y1'];
Lab2 = [z2' x2' y2'];

% Find distance between each set of consequtive points
Lab1_dif = diff(Lab1,1,1);
Lab2_dif = diff(Lab2,1,1);

% Using simple 1931 color difference metric
Lab1_inv = sqrt(sum(Lab1_dif.^2,2));
Lab2_inv = sqrt(sum(Lab2_dif.^2,2));

% Hence find length of complete arc
arclen_blu = sum(Lab1_inv);
arclen_red = sum(Lab2_inv);

% Find 256 points equally spaced out of the 100000 we generated
n_target = 256;
arclen_inv_target_blu = arclen_blu./(n_target/2-1);
arclen_inv_target_red = arclen_red./(n_target/2-1);

% Match target to the arclength NOT the euclidean length between points
% This is easier to code
Lab1_inv_runtot = [0; cumsum(Lab1_inv)];
Lab2_inv_runtot = [0; cumsum(Lab2_inv)];
arclen_targets_blu = arclen_inv_target_blu * (0:(n_target/2-1));
arclen_targets_red = arclen_inv_target_red * (0:(n_target/2-1));

t1_new = interp1(Lab1_inv_runtot, t1, arclen_targets_blu);
t2_new = interp1(Lab2_inv_runtot, t2, arclen_targets_red);

t_new = [btsp_t_srt t1_new(2:end-1) pi t2_new(2:end-1) btsp_t_end];

x = C(1) + a * cos(t_new) * U(1) + b * sin(t_new) * V(1);
y = C(2) + a * cos(t_new) * U(2) + b * sin(t_new) * V(2);
z = C(3) + a * cos(t_new) * U(3) + b * sin(t_new) * V(3);

Lab = [z' x' y'];
cform = makecform('lab2srgb');
rgb = applycform(Lab, cform);

ciebow_cmap = rgb;

img = repmat(rgb,[1 1 20]);
img = permute(img,[1 3 2]);

figure;
imagesc(img);
axis xy;
figure;
imagesc(img(1:4:end,:,:));
axis xy;

[TF,P2] = isingamut(Lab,g,'Lab');

% disp(sum(TF)/length(TF))

figure;
hold on;
plot3(Lab(:,2), Lab(:,3), Lab(:,1), 'bx')
plot3(P2(:,2) , P2(:,3) , P2(:,1) , 'r-')
xlabel('a');
ylabel('b');
zlabel('L');
view(0,90);


cmap_dif = diff(Lab,1,1);
cmap_inv = sqrt(sum(cmap_dif.^2,2));