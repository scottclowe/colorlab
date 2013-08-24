function cmap = cie_bluewhitered_cmap_old(n, verbose)

% -------------------------------------------------------------------------
% Default inputs
if nargin<1 || isempty(n)
    n = size(get(gcf,'colormap'),1);
end
if nargin<2
    verbose = false;
end

% -------------------------------------------------------------------------
% Parameters

use_uplab = false;

% CIELCH      [  L*    c    h]
  lchblue   = [  39,  83, 292];
% lchbluew  = [ 100,   0, 292];
  lchred    = [  39,  83,  40];
% lchredw   = [ 100,   0,  36];

% -------------------------------------------------------------------------

% % Finding the right values
% h = 36; %292;
% g = fetch_cielchab_gamut('srgb');
% li = g.lch(:,3)==h;
% gh = g.lch(li,:);
% figure; scatter(gh(:,2),gh(:,1),20,g.rgb(li,:));
% set(gca,'Color',[.48 .48 .48]);
% hold on;
% M = (100-gh(end-6,1))/(0-gh(end-6,2));
% C = 100; Y = 100:-1:40; X = (Y-C)/M;
% plot(X,Y,'k-');
% 
% 
% % Checking right values
% g = fetch_cielchab_gamut('srgb');
% 
% li = g.lch(:,3)==lchblue(3);
% gh = g.lch(li,:);
% figure; 
% subplot(1,2,1);
% scatter(gh(:,2),gh(:,1),20,g.rgb(li,:));
% set(gca,'Color',[.48 .48 .48]);
% hold on;
% plot([0 lchblue(2)],[100 lchblue(1)],'o-k');
% box on;
% 
% li = g.lch(:,3)==lchred(3);
% gh = g.lch(li,:);
% subplot(1,2,2);
% scatter(gh(:,2),gh(:,1),20,g.rgb(li,:));
% set(gca,'Color',[.48 .48 .48]);
% hold on;
% plot([0 lchred(2)],[100 lchred(1)],'o-k');
% box on;


% CIELab    [  L*   a*   b*]
labblue   = [lchblue(1), lchblue(2)*cos(lchblue(3)/360*2*pi), lchblue(2)*sin(lchblue(3)/360*2*pi)];
labred    = [lchred(1) ,  lchred(2)*cos(lchred(3)/360*2*pi) ,  lchred(2)*sin(lchred(3)/360*2*pi) ];

neach = ceil(n/2);

L1 = linspace(labblue(1),100,neach);
a1 = linspace(labblue(2),  0,neach);
b1 = linspace(labblue(3),  0,neach);

L2 = linspace(100,labred(1),neach);
a2 = linspace(  0,labred(2),neach);
b2 = linspace(  0,labred(3),neach);

Lab1 = [L1' a1' b1'];
Lab2 = [L2' a2' b2'];

switch mod(n,2)
    case 0
        % n is even
        Lab1 = Lab1(1:end-1,:);
        Lab2 = Lab2(2:end  ,:);
    case 1
        % n is odd
        Lab1 = Lab1(1:end-1,:);
%         Lab2 = Lab2;
end

% Convert from Lab to srgb
Lab  = [Lab1;Lab2];
cmap = soft_lab2rgb(Lab, use_uplab);


% If verbose, output figures showing colormap construction
if verbose
% Finding the right values
h = 36; %292;
g = fetch_cielchab_gamut('srgb');
li = g.lch(:,3)==h;
gh = g.lch(li,:);
figure; scatter(gh(:,2),gh(:,1),20,g.rgb(li,:));
set(gca,'Color',[.48 .48 .48]);
hold on;
M = (100-gh(end-6,1))/(0-gh(end-6,2));
C = 100; Y = 100:-1:40; X = (Y-C)/M;
plot(X,Y,'k-');


% Checking right values
g = fetch_cielchab_gamut('srgb');

% li = g.lch(:,3)==lchblue(3);
% gh = g.lch(li,:);
% figure; 
% hold on;
% subplot(1,2,1);
% % scatter(gh(:,2),gh(:,1),20,g.rgb(li,:));
% plot(gh(:,2),gh(:,1),'b-');
% set(gca,'Color',[.48 .48 .48]);
% plot([0 lchblue(2)],[100 lchblue(1)],'o-k');
% box on;
% 
% li = g.lch(:,3)==lchred(3);
% gh = g.lch(li,:);
% subplot(1,2,2);
% % scatter(gh(:,2),gh(:,1),20,g.rgb(li,:));
% plot(gh(:,2),gh(:,1),'r-');
% set(gca,'Color',[.48 .48 .48]);
% hold on;
% plot([0 lchred(2)],[100 lchred(1)],'o-k');
% box on;


figure;
hold on;
ghb = g.lch(g.lch(:,3)==lchblue(3),:);
plot(ghb(:,2),ghb(:,1),'b-');
plot([0 lchblue(2)],[100 lchblue(1)],'k-');

ghr = g.lch(g.lch(:,3)==lchred(3),:);
plot(ghr(:,2),ghr(:,1),'r-');
plot([0 lchred(2)],[100 lchred(1)],'ko');
set(gca,'Color',[.48 .48 .48]);
box on;

% Plot the colormap
img = repmat(cmap,[1 1 20]);
img = permute(img,[1 3 2]);
figure;
imagesc(img);
axis xy;
% figure;
% imagesc(img(1:4:end,:,:));
% axis xy;
end